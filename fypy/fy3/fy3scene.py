# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3scene.py

@Modify Time :  2023/10/10 17:50   

@Author      : Lee    

@Version     : 1.0   

@Description :
对FY3 MERSI L1数据进行辐射定标

针对FY-3卫星数据，包含MERSI 5分钟块、
MWRI 升降轨、MWHS、MWTS等整圈轨道数据的投影、拼接、裁剪

'''
import os
import sys
import numpy as np
import datetime
import time
import tempfile
import glob
from tqdm import tqdm
from osgeo import gdal, gdalconst, osr, ogr
from PIL import Image

from fypy.tools.hdfpro import writehdf, writehdf_fileinfo,\
                              readhdf, readhdf_fileinfo
from fypy.tools.ncpro import writenc, writenc_fileinfo
from fypy.tools.tifpro import writetiff, gettype

from fypy.fy3.fy3core import GetNameInfo, GetSourceInfo, CreateVrt,\
                             calref, calemiss, \
                             single_files, xy2latlon, latlon2xy

FY3ParmDir = os.path.relpath(os.path.dirname(__file__))

class fy3scene :

    def __init__(self, filename=None, SatID=None, InstID=None, ProdID=None,
                 RegionID=None, LevelID=None, StartTime=None, Resolution=None,
                 Proj=None, CType=None) :

        self.TempFile = []
        self.img = None

        self.Parse(filename, SatID, InstID, ProdID, RegionID,
                   LevelID, StartTime, Resolution, Proj, CType)

    def Calibration(self, filename, sdsname):
        ''' 对FY3 L1数据进行辐射定标 '''
        # EV_250_Aggr.1KM_RefSB     1~4         # 0.47, 0.55, 0.65, 0.865,
        # EV_1KM_RefSB              5~19        # 1.38, 1.64, 2.13, 0.412,
        #                                         0.443, 0.49, 0.555, 0.67,
        #                                         0.709, 0.746, 0.865, 0.905,
        #                                         0.936, 0.94, 1.03,
        # EV_1KM_Emissive           20~23       # 3.796, 4.046, 7.233, 8.56,
        # EV_250_Aggr.1KM_Emissive  24~25       # 10.714, 11.948

        data = readhdf(filename, sdsname)
        if 'EV_250_Aggr.1KM_RefSB' in sdsname :
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')
            return calref(data, coef[np.arange(1, 5)-1, :])
        elif 'EV_1KM_RefSB' in sdsname :
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')
            return calref(data, coef[np.arange(5, 20)-1, :])
        elif 'EV_1KM_Emissive' in sdsname :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return calemiss(data, 10000.0/wave_length[np.arange(20, 24)-1])
        elif 'EV_250_Aggr.1KM_Emissive' in sdsname :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return calemiss(data, 10000.0/wave_length[np.arange(24, 26)-1])
        elif sdsname in ['EV_250_RefSB_b1', 'EV_250_RefSB_b2',
                         'EV_250_RefSB_b3', 'EV_250_RefSB_b4'] :
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')
            bandid = int(sdsname[-1])
            return calref(data, [coef[bandid-1, :]])
        elif sdsname in ['EV_250_Emissive_b24', 'EV_250_Emissive_b25'] :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            bandid = int(sdsname[-2:])
            return calemiss(data, [10000.0/wave_length[bandid-1]])
        else:
            tmpdata, sdsinfo = readhdf(filename, sdsname,dictsdsinfo={})
            data = tmpdata.copy()
            if 'Slope' in sdsinfo and 'Intercept' in sdsinfo:
                data = data * sdsinfo['Slope'] + sdsinfo['Intercept']

            if 'valid_range' in sdsinfo :
                data[(tmpdata < sdsinfo['valid_range'][0]) | (tmpdata > sdsinfo['valid_range'][1])] = np.nan

            return data

    def Translate(self, srcdata, srclat, srclon, resolution=None,
                  vmin=None, vmax=None, extent=None, resampleAlg='near',
                  srcNodata=-999.0, dstNodata=None, dstSRS="EPSG:4326"):
        ''' 针对FY-3卫星数据，包含MERSI 5分钟块、MWRI 升降轨、MWHS、MWTS等整圈轨道数据的投影 '''

        if resolution is None :
            if hasattr(self, 'Resolution') :
                resolution = self.__getattribute__('Resolution')
            else:
                raise Exception('需要指定输出空间分辨率【resolution】')

        if dstNodata is None :
            dstNodata = srcNodata

        flag = (srclon > 180) | (srclon < -180) | (srclat > 90)  | (srclat < -90)
        srclon[flag] = np.nan
        srclat[flag] = np.nan
        # srcdata[flag] = np.nan

        # 获取投影数据的范围
        if extent is None :
            extent = [np.nanmin(srclon), np.nanmin(srclat),
                      np.nanmax(srclon), np.nanmax(srclat)]

        # 获取输入数据的有效值范围
        if vmin is None : vmin = np.nanmin(srcdata)
        if vmax is None : vmax = np.nanmax(srcdata)
        if np.isnan(vmin) or np.isnan(vmax) :
            print('数据有效范围为NAN，将不做投影')
            return None

        data = np.array(srcdata).copy()

        if vmax is not None and vmin is not None :
            data[(srcdata < vmin) | (srcdata > vmax)] = srcNodata

        data[np.isnan(data)] = srcNodata

        tmp_file = tempfile.NamedTemporaryFile(prefix="tmp_fypy_fy3orbit_", delete=True)
        temphdf = tmp_file.name + '.hdf'
        self.TempFile.append(temphdf)
        # 创建临时的数据文件
        writehdf(temphdf, 'data', data, overwrite=1)
        writehdf(temphdf, 'lon', srclon, overwrite=0)
        writehdf(temphdf, 'lat', srclat, overwrite=0)

        layer = GetSourceInfo(temphdf, 'data')

        vrtFile = tmp_file.name + '.vrt'
        self.TempFile.append(vrtFile)

        CreateVrt(vrtFile, temphdf, layer, '/lon', '/lat')

        ds = gdal.Warp('', vrtFile,
                         format='MEM',  geoloc=True,
                         dstSRS=dstSRS,  resampleAlg=resampleAlg,
                         srcNodata= srcNodata, dstNodata=dstNodata,
                         outputBounds=extent,  # (minX, minY, maxX, maxY)
                         xRes=resolution, yRes=resolution)

        if ds is None:
            print('处理失败')
            return None
        print('投影转换成功')

        return ds

    def Block10Merge(self, srcfile, ProdID, sdsname, fillvalue=None):
        ''' FY3极轨卫星10度块产品拼接 '''

        if isinstance(srcfile, list):
            countfile = len(srcfile)
            if countfile == 0 :
                return None
            elif countfile == 1 :
                return single_files(srcfile[0], ProdID, sdsname, fillvalue=fillvalue)
            else:
                all_ds = []

                for filename in srcfile :
                    dset = single_files(filename, ProdID, sdsname, fillvalue=fillvalue)
                    all_ds.append(dset)

                ds = gdal.Warp('', all_ds, format='MEM',)
                return ds
        elif isinstance(srcfile, str):
            return single_files(srcfile, ProdID, sdsname, fillvalue=fillvalue)

    def Hammer2Wgs84(self, srcDS, xRes=None, yRes=None, fillvalue=None, blocksize=1000):
        ''' https://blog.csdn.net/flued_g/article/details/50480508 '''

        data_src  = srcDS.ReadAsArray()
        trans     = srcDS.GetGeoTransform()
        srcNodata = srcDS.GetRasterBand(1).GetNoDataValue()
        if fillvalue is None :
            dstNodata = srcNodata
        else:
            dstNodata = fillvalue

        if xRes is None :
            xRes = trans[1]
        if yRes is None :
            yRes = trans[5]

        cols      = srcDS.RasterXSize
        rows      = srcDS.RasterYSize
        bands     = srcDS.RasterCount #波段数

        XBlockCnt = int(np.ceil((cols/blocksize)))
        YBlockCnt = int(np.ceil((rows/blocksize)))

        x_res= trans[1]
        y_res= trans[5]

        lon_max = np.nan
        lon_min = np.nan
        lat_max = np.nan
        lat_min = np.nan

        for j in np.arange(XBlockCnt) :
            for i in np.arange(YBlockCnt) :
                x_s = trans[0] + j * blocksize * x_res
                y_s = trans[3] + i * blocksize * y_res

                row = blocksize
                col = blocksize
                if i == (YBlockCnt-1) :
                    row = rows - (YBlockCnt-1) * blocksize

                if j == (XBlockCnt-1) :
                    col = cols - (XBlockCnt-1) * blocksize

                # hammer转为WGS84坐标
                lon, lat = xy2latlon(x_s, y_s, x_res, y_res, row, col)

                # 计算图像的四角范围
                lon_max =  np.nanmax([lon_max, np.nanmax(lon)])
                lon_min = np.nanmin([lon_min, np.nanmin(lon)])

                lat_max = np.nanmax([lat_max, np.nanmax(lat)])
                lat_min = np.nanmin([lat_min, np.nanmin(lat)])
                del lon, lat

                # print(lon_min, lon_max, lat_min, lat_max)

        # 计算输出图像的大小
        newcols = int(np.ceil((lon_max - lon_min) / xRes))
        newrows = int(np.ceil((lat_max - lat_min) / np.fabs(yRes)))

        # 创建输出对象
        datatype = gettype(data_src.dtype)
        driver = gdal.GetDriverByName("MEM")
        outds = driver.Create('', newcols, newrows, bands, datatype)

        # 设置仿射
        trans_dst = [lon_min, xRes, 0,
                     lat_max, 0, -yRes]
        outds.SetGeoTransform(trans_dst)

        # 设置投影坐标系
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS("EPSG:4326")
        outds.SetProjection(sr.ExportToWkt())

        XBlockCnt = int(np.ceil((newcols/blocksize)))
        YBlockCnt = int(np.ceil((newrows/blocksize)))
        with tqdm(total=XBlockCnt*YBlockCnt, iterable='iterable',
                  desc = '正在分块投影', mininterval=1) as pbar:
            for j in np.arange(XBlockCnt) :
                for i in np.arange(YBlockCnt) :

                    x_s = lon_min + j * blocksize * xRes
                    y_s = lat_max - i * blocksize * np.fabs(yRes)

                    row = blocksize
                    col = blocksize
                    if i == (YBlockCnt-1) :
                        row = newrows - (YBlockCnt-1) * blocksize

                    if j == (XBlockCnt-1) :
                        col = newcols - (XBlockCnt-1) * blocksize

                    # WGS84转为Hammer
                    x, y = latlon2xy(x_s, y_s, col, row, xRes)

                    x_index = (np.floor((x - trans[0]) / x_res)).astype(int)
                    y_index = (np.floor((y - trans[3]) / y_res)).astype(int)

                    flag = (x_index < 0) | (x_index >= cols) | \
                           (y_index < 0) | (y_index >= rows)
                    x_index[flag] = 0
                    y_index[flag] = 0
                    newdata = data_src[y_index, x_index]
                    newdata[flag] = dstNodata
                    outds.GetRasterBand(1).WriteArray(newdata, xoff=int(j*blocksize), yoff=int(i*blocksize))
                    pbar.update(1)
            pbar.close()
        if dstNodata is not None :
            outds.GetRasterBand(1).SetNoDataValue(float(dstNodata))

        return outds

    def Clip(self, srcfile, extent=None, shpname=None, resolution=None,
             dstSRS=None, epsg=None, dstNodata=None,
             resampleAlg=gdalconst.GRA_NearestNeighbour, **kwargs):
        ''' 影像裁剪 '''

        if epsg is not None and dstSRS is None:
            dstSRS = "EPSG:%d" %(int(epsg))

        if shpname is None :
            ds = gdal.Warp("", srcfile, format='MEM',
                           dstSRS=dstSRS, dstNodata=dstNodata,
                           xRes=resolution, yRes=resolution,
                           resampleAlg=resampleAlg,
                           outputBounds=extent, **kwargs)
        else:
            ds = gdal.Warp(
                "", srcfile, format='MEM', dstSRS=dstSRS,
                resampleAlg=resampleAlg, dstNodata=dstNodata,
                cutlineDSName=shpname, cropToCutline=True,
                xRes=resolution, yRes=resolution, **kwargs)

        return ds

    def load(self, filename, ProdID='truecolor'):
        data = self.Calibration(filename, '/Data/EV_250_Aggr.1KM_RefSB')
        b= data[0]
        g= data[1]
        r= data[2]

        if ProdID in ['truecolor'] :
            self.img  = self._truecolor(r, g, b)

    def show(self, ):
        if self.img is None :
            raise Exception('请先加载【load】一个对象后再【show】')

        self.img.show()

    def SaveThematic(self, outname):
        if self.img is None :
            raise Exception('请先加载【load】一个对象后再【SaveThematic】')

        self.img.save(outname)

    def DS2Tiff(self, outname, srcDS):
        ''' 保存为GeoTiff文件 '''

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()

        writetiff(outname, data, trans, prj, fillvalue=fillvalue)

    def DS2Netcdf(self, outname, sdsname, srcDS):
        ''' 保存为NetCDF文件 '''

        if srcDS is None :
            raise Exception('srcDS为None')

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()

        if len(data.shape) == 2 :
            level, (height, width) = 1, data.shape
        elif len(data.shape) == 3 :
            level, height, width = data.shape
        else:
            raise Exception('仅暂支持2D或3D数据的输出')

        lon = trans[0] + trans[1] * np.arange(width)
        lat = trans[3] + trans[5] * np.arange(height)


        dictfileinfo = {
            'trans' : trans,
            'prj'   : prj
        }
        writenc_fileinfo(outname, dictfileinfo=dictfileinfo, overwrite=1)
        writenc(outname, 'latitude',  lat, overwrite=0)
        writenc(outname, 'longitude', lon, overwrite=0)
        if level == 1 :
            writenc(outname, sdsname, data, dimension=('latitude', 'longitude'), overwrite=0)
        else:
            writenc(outname, 'level', np.arange(level), overwrite=0)
            writenc(outname, sdsname, data, dimension=('level', 'latitude', 'longitude'),
                    overwrite=0, fill_value=fillvalue)

    def DS2Hdf(self, outname, sdsname, srcDS):
        ''' 保存为HDF5文件 '''
        if srcDS is None :
            raise Exception('srcDS为None')

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()
        if data is None :
            raise Exception('读取图层数据失败')

        if len(data.shape) == 2 :
            level, (height, width) = 1, data.shape
        elif len(data.shape) == 3 :
            level, height, width = data.shape
        else:
            raise Exception('仅暂支持2D或3D数据的输出')

        lon = trans[0] + trans[1] * np.arange(width)
        lat = trans[3] + trans[5] * np.arange(height)
        dictfileinfo = {
            'trans' : trans,
            'prj'   : prj
        }
        dictsdsinfo = {
            'name' : sdsname,
            'fillvalue' : fillvalue
        }
        writenc_fileinfo(outname, dictfileinfo=dictfileinfo, overwrite=1)
        writehdf(outname,  'latitude',  lat, overwrite=0)
        writehdf(outname, 'longitude',  lon, overwrite=0)
        writehdf(outname,     sdsname, data, overwrite=0, dictsdsinfo=dictsdsinfo)

    def Parse(self, filename=None, SatID=None, InstID=None, ProdID=None,
              RegionID=None, LevelID=None, StartTime=None, Resolution=None,
              Proj=None, CType=None):
        nameinfo = {}
        if filename is not None :
            nameinfo = GetNameInfo(filename)

        if SatID is not None :
            self.SatID = SatID
        elif 'SatID' in nameinfo :
            self.SatID = nameinfo['SatID']

        if InstID is not None :
            self.InstID = InstID
        elif 'InstID' in nameinfo :
            self.InstID = nameinfo['InstID']

        if ProdID is not None :
            self.ProdID = ProdID
        elif 'ProdID' in nameinfo :
            self.ProdID = nameinfo['ProdID']

        if Proj is not None :
            self.Proj = Proj
        elif 'Proj' in nameinfo :
            self.Proj = nameinfo['Proj']

        if RegionID is not None :
            self.RegionID = RegionID
        elif 'RegionID' in nameinfo :
            self.RegionID = nameinfo['RegionID']

        if LevelID is not None :
            self.LevelID = LevelID
        elif 'LevelID' in nameinfo :
            self.LevelID = nameinfo['LevelID']

        if StartTime is not None :
            self.StartTime = StartTime
        elif 'StartTime' in nameinfo :
            self.StartTime = nameinfo['StartTime']

        if Resolution is not None :
            self.Resolution = Resolution
        elif 'Resolution' in nameinfo :
            self.Resolution = nameinfo['Resolution']

        if CType is not None :
            self.CType = CType
        elif 'CType' in nameinfo :
            self.CType = nameinfo['CType']

    def _truecolor(self, r, g, b) :

        data_b = b * 0.7
        data_g = g * 0.7
        data_r = r * 0.7

        data_b = np.array(data_b * 10000, dtype=np.int32)
        data_g = np.array(data_g * 10000, dtype=np.int32)
        data_r = np.array(data_r * 10000, dtype=np.int32)

        data_b[data_b < 0] = 0
        data_g[data_g < 0] = 0
        data_r[data_r < 0] = 0

        data_b[data_b >= 10000] = 10000 - 1
        data_g[data_g >= 10000] = 10000 - 1
        data_r[data_r >= 10000] = 10000 - 1

        if os.path.isfile(os.path.join(FY3ParmDir,'cmp.npy')):
            cmap = np.load(os.path.join(FY3ParmDir,'cmp.npy'))

        r = cmap[data_r]
        g = cmap[data_g]
        b = cmap[data_b]

        rgbArray = np.zeros((r.shape[0], r.shape[1], 4), dtype=np.float64)
        rgbArray[..., 0] = r
        rgbArray[..., 1] = g
        rgbArray[..., 2] = b
        rgbArray[..., 3] = 255

        img = Image.fromarray(rgbArray.astype(np.uint8))

        return img


    def __del__(self):

        for item in self.TempFile :
            if os.path.isfile(item) :
                try:
                    os.remove(item)
                    # print(item)
                except BaseException as e:
                    print('删除%s失败' %(item))

        # 删除系统中由fypy创建的临时文件
        tmp_file = tempfile.NamedTemporaryFile(prefix="tmp_fypy_", delete=True)
        temphdf = tmp_file.name + '.hdf'
        tempdir = os.path.dirname(temphdf)
        filelist = glob.glob(os.path.join(tempdir, 'tmp_fypy_*'))
        for filename in filelist :
            if (time.time() - os.stat(filename).st_mtime ) > 5*24*60*60 :
                if os.path.isfile(filename) :
                    try:
                        os.remove(filename)
                        # print(item)
                    except BaseException as e:
                        print('删除%s失败' %(filename))
