# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3Orbit.py

@Modify Time :  2023/8/11 17:08   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import tempfile

from osgeo import gdal, gdalconst, osr, ogr


from fypy.tools.hdfpro import writehdf
from fypy.tools.ncpro import writenc
from fypy.tools.tifpro import writetiff


class FY3Orbit() :
    '''针对FY-3卫星数据，包含MERSI 5分钟块、
    MWRI 升降轨、MWHS、MWTS等整圈轨道数据的投影、拼接、裁剪'''

    def __init__(self, satID, instID):
        self.satID = satID
        self.instID = instID

    def translate(self, srcdata, srclat, srclon, dstfile=None,
                 vmin=None, vmax=None, resolution=0.01,
                 minX=None, maxX=None, minY=None, maxY=None,  # 需要投影的区域范围
                 resampleAlg='near',
                 srcNodata=-999.0, dstNodata=None, dstSRS="EPSG:4326"):

        self.TempFile = []

        if dstNodata is None :
            dstNodata = srcNodata

        flag = (srclon > 180) | (srclon < -180) | (srclat > 90)  | (srclat < -90)
        srclon[flag] = np.nan
        srclat[flag] = np.nan
        # srcdata[flag] = np.nan

        # 获取投影数据的范围
        if maxX is None : maxX = np.nanmax(srclon)
        if minX is None : minX = np.nanmin(srclon)

        if maxY is None : maxY = np.nanmax(srclat)
        if minY is None : minY = np.nanmin(srclat)

        if vmin is None : vmin = np.nanmin(srcdata)
        if vmax is None : vmax = np.nanmax(srcdata)

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

        layer = self._GetSourceInfo(temphdf, 'data')

        vrtFile = tmp_file.name + '.vrt'
        self.TempFile.append(vrtFile)

        self._createVrt(vrtFile, temphdf, layer, '/lon', '/lat')

        if dstfile is None :
            format = 'MEM'
            dstfile = ''
            creationOptions = []
        else:
            format = 'GTiff'
            creationOptions = ["COMPRESS=LZW"]
        dset = gdal.Warp(dstfile, vrtFile,
                            format=format,  geoloc=True,
                            dstSRS=dstSRS,  resampleAlg=resampleAlg,
                            srcNodata= srcNodata, dstNodata=dstNodata,
                            outputBounds=(minX, minY, maxX, maxY),  # (minX, minY, maxX, maxY)
                            xRes=resolution, yRes=resolution,
                            creationOptions=creationOptions)

        if dset == None:
            raise Exception('处理失败')

        self.width = dset.RasterXSize #栅格矩阵的列数
        self.height = dset.RasterYSize #栅格矩阵的行数
        self.bands = dset.RasterCount #栅格矩阵的波段数

        dstdata = dset.ReadAsArray()     #获取数据
        trans = dset.GetGeoTransform()    #获取仿射矩阵信息
        prj = dset.GetProjection()       #获取投影信息

        self.dstdata = dstdata
        self.trans = trans
        self.prj = prj
        self.dstNodata = dstNodata
        self.lon = trans[0] + trans[1] * np.arange(self.width)
        self.lat = trans[3] + trans[5] * np.arange(self.height)

        print('完成投影转换')

        # if vmax is not None and vmin is not None :
        #     dstdata[(dstdata < vmin) | (dstdata > vmax)] = dstNodata

        return dstdata, trans, prj

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

    def _GetSourceInfo(self, filename, sdsname):
        src_ds = gdal.Open(filename)
        layers = src_ds.GetSubDatasets()

        # 获取sdsname所在的图层栅格索引
        src_raster = self.GetLayer(layers, sdsname)
        if src_raster is None :
            raise Exception('数据集【%s】不在文件中【%s】' %(sdsname, filename))

        return src_raster

    def GetLayer(self, layers, sdsname):
        '''
        获取指定的图层的索引名
        :param layers: tuple
        :return: str
        '''

        if sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')
                # print(self.sdsname, l_name)
                if sdsname in l_name:
                    return layer[0]

        return None

    def _createVrt(self, vrtDir, srcfile, layer, srclon=None, srclat=None):

        if layer is None :
            gdal.Translate(vrtDir,
                           srcfile,
                           format='vrt')
        else:
            gdal.Translate(vrtDir,
                           layer,
                           format='vrt')

        lines = []
        with open(vrtDir, 'r') as f:
            for line in f:
                lines.append(line)
        lines.insert(1,'<Metadata domain="GEOLOCATION">\n')
        lines.insert(2,' <MDI key="LINE_OFFSET">1</MDI>\n')
        lines.insert(3, ' <MDI key="LINE_STEP">1</MDI>\n')
        lines.insert(4, ' <MDI key="PIXEL_OFFSET">1</MDI>\n')
        lines.insert(5, ' <MDI key="PIXEL_STEP">1</MDI>\n')
        lines.insert(6, ' <MDI key="SRS">GEOGCS["WGS84",'
                        'DATUM["WGS_1984",'
                        'SPHEROID["WGS84",6378137,298.257223563,'
                        'AUTHORITY["EPSG","7030"]],'
                        'TOWGS84[0,0,0,0,0,0,0],'
                        'AUTHORITY["EPSG","6326"]],'
                        'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
                        'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],'
                        'AUTHORITY["EPSG","4326"]]</MDI>\n')
        lines.insert(7, ' <MDI key="X_BAND">1</MDI>')
        lines.insert(8, ' <MDI key="X_DATASET">HDF5:"{geofile}":/{srclon}</MDI>\n'.format(
            geofile=srcfile, srclon=srclon))
        lines.insert(9, ' <MDI key="Y_BAND">1</MDI>\n')
        lines.insert(10, ' <MDI key="Y_DATASET">HDF5:"{geofile}":/{srclat}</MDI>\n'.format(
            geofile=srcfile, srclat=srclat))
        lines.insert(11, '</Metadata>\n')
        with open(vrtDir, 'w') as f:
            for line in lines:
                f.writelines(line)

    def __del__(self):
        '''
        清理临时文件
        :return:
        '''
        for item in self.TempFile :
            if os.path.isfile(item) :
                try:
                    os.remove(item)
                    # print(item)
                except BaseException as e:
                    print('删除%s失败' %(item))


