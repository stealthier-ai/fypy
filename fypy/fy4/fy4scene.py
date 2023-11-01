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

from fypy.fy4.fy4searchtable import fy4searchtable

from fypy.tools.tifpro import writetiff
from fypy.tools.ncpro import writenc
from fypy.tools.hdfpro import writehdf


class fy4scene(fy4searchtable) :

    def __init__(self, subpoint, resolution):
        '''

        Parameters
        ----------
        subpoint: float
            星下点
        resolution: float
            数据图像的分辨率，单位为degree
        '''

        super().__init__(subpoint, resolution)

    def clip(self, data, shpname=None, outputBounds=None,
             srcNodata=65535, dstNodata=None,
             dstSRS='EPSG:4326'):
        '''
        将标称投影转换成等经纬投影（NOM->GLL）

        Parameters
        ----------
        data : numpy.array
            输入数据
        shpname: str, optional
            掩膜的面矢量（polygon）
        outputBounds : tuple, optional
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
        else:
            for i in range(im_bands):
                memDs.GetRasterBand(i+1).WriteArray(im_data[i])

        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 '
                            '+lon_0={subpoint} +no_defs'.format(
                                    subpoint=self.subpoint ))
        memDs.SetProjection(srs.ExportToWkt())
        memDs.SetGeoTransform([-5496000, self.resolution*100*1000, 0,
                               5496000, 0, -self.resolution*100*1000])

        warpDs = gdal.Warp('', memDs, format='MEM', dstSRS=dstSRS,
                           cutlineDSName=shpname, cropToCutline=True,
                           outputBounds=outputBounds, resampleAlg='near',
                           xRes=self.resolution, yRes=self.resolution,
                           srcNodata=srcNodata, dstNodata=dstNodata)

        if warpDs is None :
            return None, None, None

        dstdata = warpDs.ReadAsArray()     #获取数据
        trans = warpDs.GetGeoTransform()    #获取仿射矩阵信息
        prj = warpDs.GetProjection()       #获取投影信息

        self.width = warpDs.RasterXSize #栅格矩阵的列数
        self.height = warpDs.RasterYSize #栅格矩阵的行数
        self.bands = warpDs.RasterCount #栅格矩阵的波段数

        self.dstdata = dstdata
        self.trans = trans
        self.prj = prj
        self.dstNodata = dstNodata
        self.lon = trans[0] + trans[1] * np.arange(self.width)
        self.lat = trans[3] + trans[5] * np.arange(self.height)
        del warpDs

        return dstdata, trans, prj

    def getL1Data(self, filename, bandID=1, fillvalue=65535):
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
        Begin_Line_Number = fp.attrs['Begin Line Number'][0]
        End_Line_Number = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number = fp.attrs['End Pixel Number'][0]

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

    def getGEOData(self, filename, sdsname, fillvalue=65535):
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

    def set_green(self, vis047, vis065, fractions=(1.0, 0.13, 0.87)):
        ''' 用红光和蓝光通道模拟绿光通道 '''

        res = (vis047 * fractions[0] - vis065 * fractions[1]) / fractions[2]

        return res

    def to_tiff(self, outname):
        writetiff(outname, self.dstdata, self.trans, self.prj, fillvalue=self.dstNodata)

    def to_netcdf(self, outname):
        writenc(outname, 'latitude', self.lat, overwrite=1)
        writenc(outname, 'longitude', self.lon, overwrite=0)
        if self.bands == 1 :
            writenc(outname, 'data', self.dstdata, dimension=('latitude', 'longitude'), overwrite=0)
        else:
            writenc(outname, 'level', np.arange(self.bands), overwrite=0)
            writenc(outname, 'data', self.dstdata, dimension=('level', 'latitude', 'longitude'),
                    overwrite=0, fill_value=self.dstNodata)

    def to_hdf(self, outname):
        writehdf(outname, 'data', self.dstdata, overwrite=1)

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


