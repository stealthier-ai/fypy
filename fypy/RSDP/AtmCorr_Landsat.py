# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : AtmCorr_Landsat.py

@Modify Time :  2022/12/23 16:06   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import re
import glob
from osgeo import gdal
from Py6S import *
from .AtmCorr import AtmCorr

class AtmCorr_Landsat(AtmCorr):

    def __init__(self,  metafile):
        super(AtmCorr_Landsat, self).__init__()

        metadata = {}

        with open(metafile, 'r') as fp :
            lines = fp.readlines()

        for line in lines :
            if 'GROUP' in line :
                continue

            if '=' in line :
                line = line.replace('\n','')
                lineinfo = line.split('=')
                key = lineinfo[0].replace(' ', '')
                # value = lineinfo[1].replace(' ', '')
                value = lineinfo[1]
                metadata[key] = value

        self.metadata = metadata

    def FLAASH(self, nowdate, srcfile, BandId=None, outname=None, fillvalue=65535):

        if BandId is None :
            BandId = int((os.path.basename(srcfile).split('.')[0])[-1])

        #捕捉打开数据出错异常
        try:
            IDataSet = gdal.Open(srcfile)
        except Exception as e:
            print("文件%S打开失败" % srcfile)

        #获取行列号
        cols = IDataSet.RasterXSize
        rows = IDataSet.RasterYSize
        ImgBand = IDataSet.GetRasterBand(1)
        ImgRasterData = ImgBand.ReadAsArray(0, 0, cols, rows)

        trans = IDataSet.GetGeoTransform()
        prj = IDataSet.GetProjection()

        #辐射校正
        radiance = self.GetRadiance(BandId, ImgRasterData, self.metadata)
        #大气校正
        corrdata = self.GetAtmCorr(nowdate, BandId, radiance, self.metadata, dem=0)

        print('第%s波段大气校正完成' %BandId)

        return corrdata, trans, prj


    def GetAtmCorr(self, nowdate, BandId, radiance, metadata, dem):
        '''

        :param nowdate:
        :param BandId:
        :param radiance:
        :param metadata:
        :param dem:
        :return:
        '''
        # minwl, maxwl, response = atmc.getresp(respfile)

        # 中心经纬度
        point1lat = float(metadata['CORNER_UL_LAT_PRODUCT'])
        point1lon = float(metadata['CORNER_UL_LON_PRODUCT'])
        point2lat = float(metadata['CORNER_UR_LAT_PRODUCT'])
        point2lon = float(metadata['CORNER_UR_LON_PRODUCT'])
        point3lat = float(metadata['CORNER_LL_LAT_PRODUCT'])
        point3lon = float(metadata['CORNER_LL_LON_PRODUCT'])
        point4lat = float(metadata['CORNER_LR_LAT_PRODUCT'])
        point4lon = float(metadata['CORNER_LR_LON_PRODUCT'])

        sLongitude = (point1lon + point2lon + point3lon + point4lon) / 4
        sLatitude = (point1lat + point2lat + point3lat + point4lat) / 4

        sunz = 90-float(metadata['SUN_ELEVATION'])
        suna =    float(metadata['SUN_AZIMUTH'])

        # 通过研究去区的范围去求DEM高度。
        # pointUL = dict()
        # pointDR = dict()
        # pointUL["lat"] = point1lat
        # pointUL["lon"] = point1lon
        # pointDR["lat"] = point4lat
        # pointDR["lon"] = point2lon

        # 校正波段（根据波段名称）
        if BandId == 1:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B1

        elif BandId == 2:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B2

        elif BandId == 3:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B3

        elif BandId == 4:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B4

        elif BandId == 5:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B5

        elif BandId == 6:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B6

        elif BandId == 7:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B7

        elif BandId == 8:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B8

        elif BandId == 9:
            wavelength = PredefinedWavelengths.LANDSAT_OLI_B9

        self.set_geom(nowdate, sunz=sunz, suna=suna, satz=0, sata=0)
        self.set_atm(sLatitude=sLatitude, nowdate=nowdate)
        self.set_aer()
        self.set_vis()
        self.set_altitude(dem=dem)
        self.set_resp(wavelength=wavelength)
        a, b, c  = self.corrCoeff()
        corrdata = self.corrImage(radiance, a, b, c)

        return np.array(corrdata)


    # 逐波段辐射定标
    def GetRadiance(self, BandId, data, metadata, fillvalue=65535):
        ''' 计算辐射亮度参数：DN -> radiance '''

        Gain = float(metadata['RADIANCE_MULT_BAND_%d' %(BandId)])
        Bias = float(metadata['RADIANCE_ADD_BAND_%d' %(BandId)])

        RaCal = np.where(data>0 ,Gain * data + Bias, fillvalue)

        return RaCal

    def GetReflectance(self, BandId, data, metadata, fillvalue=65535):
        ''' 计算反射率：DN -> Reflectance '''

        Gain = float(metadata['REFLECTANCE_MULT_BAND_%d' %(BandId)])
        Bias = float(metadata['REFLECTANCE_ADD_BAND_%d' %(BandId)])

        RaCal = np.where(data>0 ,Gain * data + Bias, fillvalue)

        return RaCal
