# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : AtmCorr_Sentinel.py

@Modify Time :  2023/4/27 14:47   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import json
from tqdm import tqdm
from osgeo import gdal, osr, ogr
from dateutil.relativedelta import relativedelta
import xml.dom.minidom
from Py6S import *

from .AtmCorr import AtmCorr, reflectance2radiance

EXEPATH = os.path.split(os.path.realpath(__file__))[0]

class AtmCorr_Sentinel(AtmCorr):

    def __init__(self):
        super(AtmCorr_Sentinel, self).__init__()

    def FLAASH(self, outname, tiffFile, metadata, SatID, InstID, BandId,
               fillvalue=65535, blocksize=1000):

        outdir = os.path.dirname(outname)
        if not os.path.isdir(outdir) :
            os.makedirs(outdir)
            print('成功创建路径【%s】' %(outdir))

        if not os.path.isfile(tiffFile) :
            print('文件不存在【%s】' %(tiffFile))
            return None

        # 捕捉打开数据出错异常
        try:
            srcdataset = gdal.Open(tiffFile)
        except Exception as e:
            print("文件{file}打开失败".format(file=tiffFile))

        if srcdataset == None:
            print("{file}数据集读取为空".format(file=tiffFile))
            return None

        # 获取行列号
        cols = srcdataset.RasterXSize                  # 栅格矩阵的列数
        rows = srcdataset.RasterYSize                  # 栅格矩阵的行数
        Bands = srcdataset.RasterCount                 # 波段数

        srcTrans = srcdataset.GetGeoTransform()   # 获取仿射矩阵信息
        srcProj = srcdataset.GetProjection()             # 获取投影信息

        # 创建输出结果文件
        driver = gdal.GetDriverByName("GTiff")
        dstdataset = driver.Create(outname, cols, rows, Bands, gdal.GDT_UInt16,
                                   options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdataset.SetGeoTransform(srcTrans)
        dstdataset.SetProjection(srcProj)

        # 表观反射率转换为辐射亮度值
        CalCoeffFile = os.path.join(EXEPATH, "resp", "Sentinel",
                                    "CalibrationCoefficient.json")
        if not os.path.isfile(CalCoeffFile) :
            raise Exception('哨兵卫星定标系数文件不存在【%s】' %(CalCoeffFile))
        CalCoeff = json.load(open(CalCoeffFile))

        satid = SatID.replace('-', '')
        instid = InstID.replace('-', '')

        if not satid in CalCoeff :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(SatID), CalCoeff.keys())

        if not instid in CalCoeff[satid] :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(InstID), CalCoeff[satid].keys())

        Esun = CalCoeff[SatID][InstID]['ESUN'][BandId-1]

        #大气校正
        self.setParam(metadata['Date'], metadata, SatID, InstID, BandId, dem=0)
        Coeffa, Coeffb, Coeffc  = self.corrCoeff()

        # 创建输出结果文件
        driver = gdal.GetDriverByName("GTiff")
        dstdataset = driver.Create(outname, cols, rows, Bands, gdal.GDT_UInt16,
                                   options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdataset.SetGeoTransform(srcTrans)
        dstdataset.SetProjection(srcProj)

        ReadBand = srcdataset.GetRasterBand(1)
        outband = dstdataset.GetRasterBand(1)
        outband.SetNoDataValue(fillvalue)

        #进度条参数
        XBlockcount = np.ceil(cols / blocksize)
        YBlockcount = np.ceil(rows / blocksize)
        i = 0
        j = 0
        try:
            with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                      desc = '正在进行第%i波段校正' %(BandId), mininterval=1) as pbar:
                while i < rows:
                    while j < cols:
                        # 保存分块大小
                        nXBK = blocksize
                        nYBK = blocksize

                        # 最后一块
                        if i+blocksize > rows:
                            nYBK = rows - i
                        if j+blocksize > cols:
                            nXBK=cols - j

                        # 分块读取影像
                        srcdata = ReadBand.ReadAsArray(j, i, nXBK,nYBK)
                        TOARef = srcdata * 0.0001
                        radiance = reflectance2radiance(TOARef, Esun, metadata['sunz'], metadata['Date'])

                        outImage =np.where(radiance>0, radiance, fillvalue)

                        y = np.where(outImage!=fillvalue, Coeffa * outImage - Coeffb, fillvalue)
                        atcImage = np.where(y!=fillvalue, (y / (1 + y * Coeffc))*10000, fillvalue)

                        outband.WriteArray(atcImage, j, i)
                        j=j+nXBK
                        pbar.update(1)
                    j=0
                    i=i+nYBK
        except BaseException :
            pbar.close()

        pbar.close()

    def setParam(self, nowdate, metadata, SatID, InstID, BandId, dem=0.010):
        '''

        :param nowdate:
        :param BandId:
        :param radiance:
        :param metadata:
        :param dem:
        :return:
        '''

        self.set_geom(nowdate, sunz=metadata['sunz'], suna=metadata['suna'], satz=0, sata=0)
        self.set_atm(sLatitude=metadata['centre_pos'][1], nowdate=nowdate)
        self.set_aer()
        self.set_vis()
        self.set_altitude(dem=dem)

        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        SRFFile = os.path.join(EXEPATH, 'resp', 'Sentinel', "Sentinel.json")
        if not os.path.isfile(SRFFile) :
            raise Exception('哨兵卫星光谱响应文件不存在【%s】' %(SRFFile))
        SRF = json.load(open(SRFFile))

        minwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][0]
        maxwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][1]
        response = SRF[SatID][InstID]['B%d' %(BandId)]['SRF']
        self.set_resp(minwl, maxwl, response)

        # self.set_resp(wavelength=wavelength)

    def GetMeta(self, metafile):
        dom = xml.dom.minidom.parse(metafile)
        dict_meta = {}

        #太阳天顶角、方位角
        SunAngle = dom.getElementsByTagName('Mean_Sun_Angle')
        dict_meta['sunz'] = float(SunAngle[0].getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
        dict_meta['suna'] = float(SunAngle[0].getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

        #卫星天顶角、方位角
        ViewAngles = dom.getElementsByTagName('Mean_Viewing_Incidence_Angle')
        dict_satz = {}
        dict_sata = {}

        for angle in ViewAngles:
            ViewAngle = int(angle.getAttribute('bandId'))
            dict_satz[ViewAngle+1] = float(angle.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
            dict_sata[ViewAngle+1]= float(angle.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

        dict_meta["ViewZenithAngle"] = dict_satz
        dict_meta["ViewAzimuthAngle"] = dict_sata
        # 日期:月、日  2023-04-11T03:50:21.042436Z
        dict_meta["Date"] = datetime.datetime.strptime(dom.getElementsByTagName('SENSING_TIME')[0].firstChild.data,
                                                       '%Y-%m-%dT%H:%M:%S.%fZ')

        #求影像中心经纬度
        PointULX = int(dom.getElementsByTagName('ULX')[0].firstChild.data)
        PointULY = int(dom.getElementsByTagName('ULY')[0].firstChild.data)

        Imgsizes = dom.getElementsByTagName('Size')

        for Imgsize in Imgsizes:
            Resolution = Imgsize.getAttribute('resolution')
            if Resolution == '10':
                dict_meta["Nrows"] = int(Imgsize.getElementsByTagName('NROWS')[0].firstChild.data)
                dict_meta["Ncols"] = int(Imgsize.getElementsByTagName('NCOLS')[0].firstChild.data)

        PointBRX = PointULX + 10*dict_meta["Ncols"]
        PointBRY = PointULY - 10*dict_meta["Nrows"]

        # 将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        Proj = dom.getElementsByTagName('HORIZONTAL_CS_CODE')[0].firstChild.data
        ProjCode = int(Proj.split(':')[1])

        source = osr.SpatialReference()
        source.ImportFromEPSG(ProjCode)
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        ct = osr.CoordinateTransformation(source,target)
        CoordsUL,CoordsBR = ct.TransformPoints([(PointULX,PointULY),(PointBRX,PointBRY)])

        ULLat = CoordsUL[0]
        ULLon = CoordsUL[1]
        BRLat = CoordsBR[0]
        BRLon = CoordsBR[1]

        sLongitude = (ULLon+BRLon) / 2
        sLatitude = (ULLat+BRLat) / 2

        dict_meta['extent'] = [ULLon, BRLat, BRLon, ULLat]  # minX, minY, maxX, maxY
        dict_meta['centre_pos'] = [sLongitude, sLatitude]

        return dict_meta