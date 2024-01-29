# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4pro.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee

@Version     : 1.0

@Description :

'''
import os
import sys
import numpy as np
from osgeo import gdal, osr

from fypy.fy4.fy4core import fy4searchtable

from fypy.tools.tifpro import writetiff
from fypy.tools.ncpro import writenc, writenc_fileinfo
from fypy.tools.hdfpro import writehdf

from fypy.fy4.fy4core import GetNameInfo
from fypy.fy4.fy4config import AreaInfo
from PIL import Image
RATE = 4


class fy4scene(fy4searchtable) :

    def __init__(self, filename=None, satID=None, instID=None, prodID=None,
                 sublon=None, resolution=None, regionID=None, levelID=None,
                 startTime=None, endTime=None, proj=None):
        ''' 通过文件名获取文件相关信息 '''

        self.Parse(filename, satID, instID, prodID, sublon, regionID,
                   levelID, startTime, endTime, resolution, proj)

        if (not hasattr(self, 'SubLon')) or (not hasattr(self, 'Resolution')) :
            raise Exception('请设置【sublon】和【resolution】参数')

        super().__init__(self.SubLon, self.Resolution)

    def Clip(self, data, shpname=None, extent=None,
             srcNodata=65535, dstNodata=None, dstSRS='EPSG:4326',
             resampleAlg='near'):
        '''
        将标称投影转换成等经纬投影（NOM->GLL）

        Parameters
        ----------
        data : numpy.array
            输入数据
        shpname: str, optional
            掩膜的面矢量（polygon）
        extent : list, optional
            output bounds as (minX, minY, maxX, maxY) in target SRS
        dstNodata : float
            数据填充值
        dstSRS : str
            输出投影坐标系

        Returns
        -------

        '''

        im_data = np.array(data, dtype=np.float32)
        dtype = self._gettype(im_data.dtype)

        # 只支持2、3维数据处理
        if len(im_data.shape) == 2:
            im_bands, (im_height, im_width) = 1,im_data.shape
        elif len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape

        Driver = gdal.GetDriverByName('MEM')
        memDs = Driver.Create('', im_height, im_width, im_bands, dtype)

        # 写入数据
        if im_bands == 1:
            memDs.GetRasterBand(1).WriteArray(im_data)
            if srcNodata is not None:
                memDs.GetRasterBand(1).SetNoDataValue(float(srcNodata))
        else:
            for i in range(im_bands):
                memDs.GetRasterBand(i+1).WriteArray(im_data[i])
                if srcNodata is not None:
                    memDs.GetRasterBand(i+1).SetNoDataValue(float(srcNodata))

        # 设置参考投影
        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 '
                            '+lon_0={sublon} +no_defs'.format(sublon=self.SubLon))
        memDs.SetProjection(srs.ExportToWkt())

        # 设置仿射变换
        memDs.SetGeoTransform([-5496000, self.Resolution*100*1000, 0,
                               5496000, 0, -self.Resolution*100*1000])

        # 影像裁剪
        warpDs = gdal.Warp('', memDs, format='MEM', dstSRS=dstSRS,
                           cutlineDSName=shpname, cropToCutline=True,
                           outputBounds=extent, resampleAlg=resampleAlg,
                           xRes=self.Resolution, yRes=self.Resolution,
                           srcNodata=srcNodata, dstNodata=dstNodata)

        if warpDs is None :
            return None

        return warpDs

    def Calibration(self, filename, bandID=1, fillvalue=65535):
        '''
        读取FY4 L1数据，并完成辐射定标
        Parameters
        ----------
        filename: str
            L1数据文件名
        bandID : int
            波段索引，从1开始

        Returns
        -------
            numpy.array
            辐射定标转换后的ref或bt
        '''
        if not os.path.isfile(filename) :
            print('文件不存在【%s】' %(filename))
            return None

        import h5py
        fp = h5py.File(filename, 'r')
        # 转换到区域的行列号（考虑去除图像偏移）
        Begin_Line_Number  = fp.attrs['Begin Line Number'][0]
        End_Line_Number    = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number   = fp.attrs['End Pixel Number'][0]

        cal = fp['CALChannel%02d' %(bandID)][:]
        dn = fp['NOMChannel%02d' %(bandID)][:]
        fp.close()

        flag = dn>=len(cal)
        dn[flag] = 0

        data1 = cal[dn]
        data1[flag] = fillvalue

        data = np.full(shape=(self.rowmax, self.colmax),
                       dtype=np.float32,
                       fill_value=fillvalue)
        data[Begin_Line_Number:End_Line_Number+1, Begin_Pixel_Number:End_Pixel_Number+1] = data1

        return data

    def GetGEOData(self, filename, sdsname, fillvalue=65535):
        import h5py
        fp = h5py.File(filename, 'r')
        data1 = fp[sdsname][:]

        # 转换到区域的行列号（考虑去除图像偏移）
        Begin_Line_Number = fp.attrs['Begin Line Number'][0]
        End_Line_Number = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number = fp.attrs['End Pixel Number'][0]
        fp.close()

        data = np.full(shape=(self.rowmax, self.colmax),
                       dtype=np.float32,
                       fill_value=fillvalue)
        data[Begin_Line_Number:End_Line_Number+1, Begin_Pixel_Number:End_Pixel_Number+1] = data1

        return data

    def load(self, filename, ProdID=None):
        vis045 = self.Calibration(filename, bandID=1)
        vis065 = self.Calibration(filename, bandID=2)
        vis085 = self.Calibration(filename, bandID=3)


        vis055 = self.set_green(vis045, vis065)

        rr = np.array(vis065*255, dtype=np.uint8)
        gg = np.array(vis055*255, dtype=np.uint8)
        bb = np.array(vis045*255, dtype=np.uint8)

        vm = Image.merge('RGB', [self.Arr2Img(i) for i in (rr, gg, bb)])
        vm.save(r'D:\DATA\FY4A\test.png', quality=90)

    def Arr2Img(self, arr):
        if arr.dtype == np.uint8:
            return Image.fromarray(arr, "L")
        if arr.dtype == np.uint16:
            return Image.fromarray(self.T0_255(arr), "L")
        return Image.fromarray(self.T0_255(arr.astype('u2')), "L")

    def T0_255(self, raw):
        return (raw >> RATE).astype(np.uint8)

    def show(self, filename, ProdID=None):
        pass

    def SaveThematic(self, outname, srcfile, ProdID=None):
        pass

    def set_green(self, vis047, vis065, fractions=(1.0, 0.13, 0.87)):
        ''' 用红光和蓝光通道模拟绿光通道 '''

        res = (vis047 * fractions[0] - vis065 * fractions[1]) / fractions[2]

        return res

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

    def Parse(self, filename=None, SatID=None, InstID=None, ProdID=None, SubLon=None,
              RegionID=None, LevelID=None, StartTime=None, EndTime=None, Resolution=None,
              Proj=None):
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

        if EndTime is not None :
            self.EndTime = EndTime
        elif 'EndTime' in nameinfo :
            self.EndTime = nameinfo['EndTime']

        if Resolution is not None :
            self.Resolution = Resolution
        elif 'Resolution' in nameinfo :
            self.Resolution = nameinfo['Resolution']

        if SubLon is not None :
            self.SubLon = SubLon
        elif 'SubLon' in nameinfo :
            self.SubLon = float(nameinfo['SubLon'].replace('E', ''))/10.0

    def _gettype(self, datatype):
        ''' 根据numpy的数据类型，匹配GDAL中的数据类型 '''

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
        elif datatype == np.float32 or datatype.str in ['>f4', '<f4']:
            return gdal.GDT_Float32
        elif datatype == np.float64 or datatype.str in ['>f8', '<f8']:
            return gdal.GDT_Float64
        else:
            return gdal.GDT_Unknown


