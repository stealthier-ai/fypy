# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3Block10.py

@Modify Time :  2023/8/11 17:07   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
from tqdm import tqdm


from osgeo import gdal, gdalconst, osr, ogr

from fypy.tools.hdfpro import readhdf
from fypy.tools.tifpro import writetiff, getDateSet, gettype

from .fy3config import FY3Block10CoefX, FY3Block10CoefY, FY3ProdInfo, Prj_Info



class FY3Block10():

    def __init__(self, satID=None, instID=None):
        self.SatID = satID
        self.InstID = instID

    def merge_block(self, srcfile, productName, sdsname):
        if isinstance(srcfile, list):
            self.dset = self.multiple_files(srcfile, productName, sdsname)
        elif isinstance(srcfile, str):
            self.dset = self.single_files(srcfile, productName, sdsname)

    def get_nameinfo(self, filename):
        ''' 根据输入数据的文件名获取信息 '''

        dict_fileinfo = {}

        basename = os.path.basename(filename)
        if len(basename) == 57 :
            basename = basename.split('.')[0]
            names = basename.split('_')
            dict_fileinfo['satID'] = names[0]
            dict_fileinfo['instID'] = names[1]
            dict_fileinfo['blockID'] = names[2]
            dict_fileinfo['levelID'] = names[3]
            dict_fileinfo['prodID'] = names[4]
            dict_fileinfo['prjType'] = names[6]
            dict_fileinfo['datetime'] = names[7]
            dict_fileinfo['peroid'] = names[8]
            dict_fileinfo['peroid'] = names[9]
            if 'CM' in names[9] :
                Res = int(names[9].replace('CM','')) / 10000000
            elif 'KM' in names[9] :
                Res = int(names[9].replace('KM','')) / 100
            else:
                Res = int(names[9].replace('M','')) / 100000

            dict_fileinfo['resolution'] = Res
            dict_fileinfo['strRes'] = names[9]
        else:
            print('输入文件名为非标准文件名，将不作信息提取【%s】' %(filename))
            return None

        return dict_fileinfo

    def single_files(self, filename, productName, sdsname) :
        '''单个数据产品投影拼接'''

        NameInfo = self.get_nameinfo(filename)
        if NameInfo is None :
            raise Exception('请传入官方标准文件名格式')

        resolution = NameInfo['resolution']

        if self.SatID is None :
            self.SatID = NameInfo['satID']

        if self.InstID is None :
            self.InstID = NameInfo['instID']

        prjType = NameInfo['prjType']
        if prjType in Prj_Info :
            prjwkt = Prj_Info[prjType]
        else:
            prjwkt = Prj_Info['GLL']

        blockID = NameInfo['blockID']

        lat = FY3Block10CoefY[blockID[0:2]]
        lon = FY3Block10CoefX[blockID[2:4]]

        mtrans = (float(lon) * 100000, float(resolution) * 100000, 0,
                  float(lat) * 100000, 0, -1 * float(resolution) * 100000)

        data = self.get_data(filename, productName, sdsname)
        ds = getDateSet(data, mtrans, prj=prjwkt, fillvalue=self.fillvalue)

        return ds


    def xy2latlon(self, x_s, y_s, x_res, y_res, rows, cols):
        ''' hammer转为WGS84坐标 '''
        # z^2 = 1 - x^2/2 - y^2/2
        # longitude = 2 * atan(sqrt(2) * x * z / (2 * z^2 - 1))
        # latitude = asin(sqrt(2) * y * z)

        x, y = np.meshgrid(np.linspace(x_s, x_s+x_res*(cols-1), num=cols),
                           np.linspace(y_s, y_s+y_res*(rows-1), num=rows))

        # 将degree转为meter: 100.0 * 1000.0
        x = np.where(x > (180.0 * 100.0 * 1000.0), (180.0 * 100.0 * 1000.0) - x, x)
        x = x / (180.0 * 100.0 * 1000.0)

        y = np.where(y > (90.0 * 100.0 * 1000.0), (90.0 * 100.0 * 1000.0) - y, y)
        y = y / (9000.0 * 1000.0)

        z = np.sqrt(1 - np.square(x) / 2.0 - np.square(y) / 2.0)
        lon = 2 * np.arctan(np.sqrt(2) * x * z / (2.0 * (np.square(z)) - 1))
        lat = np.arcsin(np.sqrt(2) * y * z)

        # 弧度转为度
        lon = lon / np.pi * 180.0
        lat = lat / np.pi * 180.0

        return lon, lat

    def latlon2xy(self, x, y, Pixel, Line, resolution):

        lon, lat = np.meshgrid(np.linspace(x, x+resolution*(Pixel-1), num=Pixel),
                                         np.linspace(y,y-1 * resolution*(Line-1), num=Line))

        # 将经纬度坐标转化为Hammer坐标
        lon = np.where(lon > 180.0, lon - 360.0, lon)
        lon = lon / 180.0 * np.pi
        lat = lat / 180.0 * np.pi

        newz = np.sqrt(1 + np.cos(lat) * np.cos(lon / 2.0))
        x = np.cos(lat) * np.sin(lon / 2.0) / newz
        y = np.sin(lat) / newz

        del lat, lon

        x = x * (180.0 * 100.0 * 1000.0)
        y = y * (90.0  * 100.0 * 1000.0)

        return x, y

    def hammer2wgs84(self, ds, resolution, blocksize=1000):
        ''' https://blog.csdn.net/flued_g/article/details/50480508 '''

        data_src = ds.ReadAsArray()
        trans = ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        bands = ds.RasterCount #波段数

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
                lon, lat = self.xy2latlon(x_s, y_s, x_res, y_res, row, col)

                # 计算图像的四角范围
                lon_max =  np.nanmax([lon_max, np.nanmax(lon)])
                lon_min = np.nanmin([lon_min, np.nanmin(lon)])

                lat_max = np.nanmax([lat_max, np.nanmax(lat)])
                lat_min = np.nanmin([lat_min, np.nanmin(lat)])
                del lon, lat

                # print(lon_min, lon_max, lat_min, lat_max)

        # 计算输出图像的大小
        newcols = int(np.ceil((lon_max - lon_min) / resolution))
        newrows = int(np.ceil((lat_max - lat_min) / resolution))

        # 创建输出对象
        datatype = self.gettype(data_src.dtype)
        driver = gdal.GetDriverByName("MEM")
        outds = driver.Create('', newcols, newrows, bands, datatype)

        # 设置仿射
        trans_dst = [lon_min, resolution, 0,
                     lat_max, 0, -resolution]
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

                    x_s = lon_min + j * blocksize * resolution
                    y_s = lat_max - i * blocksize * resolution

                    row = blocksize
                    col = blocksize
                    if i == (YBlockCnt-1) :
                        row = newrows - (YBlockCnt-1) * blocksize

                    if j == (XBlockCnt-1) :
                        col = newcols - (XBlockCnt-1) * blocksize

                    # WGS84转为Hammer
                    x, y = self.latlon2xy(x_s, y_s, col, row, resolution)

                    x_index = (np.floor((x - trans[0]) / x_res)).astype(int)
                    y_index = (np.floor((y - trans[3]) / y_res)).astype(int)

                    flag = (x_index < 0) | (x_index >= cols) | \
                           (y_index < 0) | (y_index >= rows)
                    x_index[flag] = 0
                    y_index[flag] = 0
                    newdata = data_src[y_index, x_index]
                    newdata[flag] = self.fillvalue
                    outds.GetRasterBand(1).WriteArray(newdata, xoff=int(j*blocksize), yoff=int(i*blocksize))
                    pbar.update(1)
            pbar.close()
        outds.GetRasterBand(1).SetNoDataValue(float(self.fillvalue))

        data = outds.ReadAsArray()
        trans = outds.GetGeoTransform()
        prj = outds.GetProjection()

        return data, trans, prj

    def gettype(self, datatype):
        '''
        根据numpy的数据类型，匹配GDAL中的数据类型
        :param datatype:
        :return: GDAL数据类型
        '''

        if datatype == np.byte or datatype == np.uint8:
            return gdal.GDT_Byte
        elif datatype == np.uint16 :
            return gdal.GDT_UInt16
        elif datatype == np.int16 :
            return gdal.GDT_Int16
        elif datatype == np.uint32 :
            return gdal.GDT_UInt32
        elif datatype == np.int32 :
            return gdal.GDT_Int32
        elif datatype == np.float32 :
            return gdal.GDT_Float32
        elif datatype == np.float64 :
            return gdal.GDT_Float64
        else:
            return gdal.GDT_Unknown

    def save_tiff(self, outname):

        data = self.dset.ReadAsArray()
        trans = self.dset.GetGeoTransform()
        prj = self.dset.GetProjection()

        writetiff(outname=outname, data=data, trans=trans, prj=prj, fillvalue=self.fillvalue)

    def clip(self, extent=None, shpname=None, resolution=None,
             dstSRS=None, epsg=None, dstNodata=None,
             resampleAlg=gdalconst.GRA_NearestNeighbour, **kwargs):
        ''' 影像裁剪 '''

        if epsg is not None and dstSRS is None:
            dstSRS = "EPSG:%d" %(int(epsg))

        if shpname is None :
            ds = gdal.Warp("", self.dset, format='MEM',
                           dstSRS=dstSRS, dstNodata=dstNodata,
                           xRes=resolution, yRes=resolution,
                           resampleAlg=resampleAlg,
                           outputBounds=extent, **kwargs)
        else:
            ds = gdal.Warp(
                "", self.dset, format='MEM', dstSRS=dstSRS,
                resampleAlg=resampleAlg, dstNodata=dstNodata,
                cutlineDSName=shpname, cropToCutline=True,
                xRes=resolution, yRes=resolution, **kwargs)

        data = ds.ReadAsArray()
        trans = ds.GetGeoTransform()      # 获取仿射矩阵信息
        prj = ds.GetProjection()          # 获取投影信息

        im_bands = ds.RasterCount #波段数
        # 对填充值进行NAN填充
        for iband in range(im_bands) :
            fillvalue = ds.GetRasterBand(iband+1).GetNoDataValue()
            if data.dtype == np.float16 or data.dtype == np.float32 or data.dtype == np.float64 :
                data[data==fillvalue] = np.nan

        ds = None
        return data, trans, prj

    def multiple_files(self, filelist, productName, sdsname) :
        '''多个数据产品投影拼接'''

        countfile = len(filelist)

        if countfile == 0 :
            return None
        elif countfile == 1 :
            self.single_files(filelist[0], productName, sdsname)
        else:
            all_ds = []

            for filename in filelist :
                dset = self.single_files(filename, productName, sdsname)
                all_ds.append(dset)

            ds = gdal.Warp('', all_ds, format='MEM',)
                           # xRes=self.resolution, yRes=self.resolution,
                           # srcNodata=self.fillvalue, dstNodata=self.fillvalue)

            return ds

    def get_bounding_box(self, filelist):

        tileindex = []

        for filename in filelist :
            basename = os.path.basename(filename)
            names = basename.split('_')
            blockID = names[2]

            lat = FY3Block10CoefY[blockID[0:2]]
            lon = FY3Block10CoefX[blockID[2:4]]

            tileindex.append([lat, lon])

        tileindex = np.array(tileindex)
        maxY = np.nanmax(tileindex[:,0])
        minY = np.nanmin(tileindex[:,0])

        maxX = np.nanmax(tileindex[:,1])
        minX = np.nanmin(tileindex[:,1])

        extent = [minX, maxX+10, minY, maxY+10]

        return extent

    def get_data(self, filename, productName, sdsname):

        try:
            if productName in FY3ProdInfo[self.SatID][self.InstID]:
                val, sdsinfo = readhdf(filename, sdsname, dictsdsinfo={})
                if 'FillValue' in sdsinfo :

                    if isinstance(sdsinfo['FillValue'], np.ndarray) :
                        self.fillvalue = sdsinfo['FillValue'][0]
                    else:
                        self.fillvalue = sdsinfo['FillValue']
                else:
                    self.fillvalue = 65535

                if 'valid_range' in sdsinfo :
                    vmin = sdsinfo['valid_range'][0]
                    vmax = sdsinfo['valid_range'][1]
                    fillflag = (val<vmin) | (val>vmax)

                if 'Slope' in sdsinfo and 'Intercept' in sdsinfo :
                    val = val * sdsinfo['Slope'] + sdsinfo['Intercept']

                val[fillflag] = self.fillvalue
                # if 'valid_range' in sdsinfo :
                #     vmin = sdsinfo['valid_range'][0]
                #     vmax = sdsinfo['valid_range'][1]
                    # fillflag = (val<vmin) | (val>vmax)
                    # val[fillflag] = self.fillvalue

                if productName == 'CLM':
                    data_u8 = np.uint8(val)
                    data_bit = np.unpackbits(data_u8).reshape(data_u8.size, 8)
                    bit21 = data_bit[:, 5] * 4 + data_bit[:, 6] * 2 + data_bit[:, 7]
                    val = bit21.reshape(val.shape)
                    # 001 = Cloudy （1）
                    # 011 = Uncertain （3）
                    # 101 = Probably  Clear （5）
                    # 111 = Confident  Clear （7）
                    val[(val != 1) & (val != 3) & (val != 5) & (val != 7)] = 255
                    val[val==1] = 0
                    val[val==3] = 1
                    val[val==5] = 2
                    val[val==7] = 3
                    self.fillvalue = 255

                key = sdsname.replace(' ', '_')
                if '/' in key :
                    key = key.split('/')[-1]

                return val
            else:
                print('产品【%s】不在绘图要求列表之内' %(productName))
                return None
        except BaseException as  e :
            print('读取【%s】：%s失败！！！' %(productName, filename))
            return None

    def __del__(self):

        self.dset = None

