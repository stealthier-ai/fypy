# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : AtmCorr_FY3D_MERSI.py

@Modify Time :  2022/12/23 15:46   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import json
from .AtmCorr import AtmCorr, reflectance2radiance
from fypy.tools import readhdf, readhdf_fileinfo, writehdf
from scipy import interpolate


EXEPATH = os.path.split(os.path.realpath(__file__))[0]


class AtmCorr_FY3D_MERSI(AtmCorr):

    def __init__(self):

        super(AtmCorr_FY3D_MERSI, self).__init__()


    def FLAASH(self, nowdate, srcdata, lat, lon,
               suna, sunz, sata, satz,
               SatID, InstID, BandId,
               dem=0.001, fillvalue=65535):
        '''
        基于Py6S进行大气校正

        Parameters
        ----------
        nowdate ： datetime
            数据时间
        srcdata : array
            输入数据
        lat : array
            纬度
        lon : array
            经度
        suna ：
            太阳方位角，degree
        sunz ：
            太阳天顶角，degree
        sata ：
            卫星方位角， degree
        satz ：
            卫星天顶角， degree
        dem ：
            高程， 单位：米
        SatID ：str
            卫星ID， FY3D
        InstID ：str
            仪器ID， MERSI
        BandId ：int
            波段ID，从1开始
        fillvalue ：float
            填充值

        Returns
        -------
            订正后的反射率
        '''
        # 表观反射率转换为辐射亮度值
        CalCoeffFile = os.path.join(EXEPATH, "resp", "FY",
                                    "CalibrationCoefficient.json")
        if not os.path.isfile(CalCoeffFile) :
            raise Exception('风云卫星定标系数文件不存在【%s】' %(CalCoeffFile))
        CalCoeff = json.load(open(CalCoeffFile))

        satid = SatID.replace('-', '')
        instid = InstID.replace('-', '')

        if not satid in CalCoeff :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(SatID), CalCoeff.keys())

        if not instid in CalCoeff[satid] :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(InstID), CalCoeff[satid].keys())


        Esun = CalCoeff[SatID][InstID]['ESUN'][BandId-1]
        # MW = CalCoeff[SatID][InstID]['MW'][BandId-1] * 0.001
        radiance = reflectance2radiance(srcdata, Esun, sunz, nowdate)

        radiance = radiance*10
        # 模型参数设置
        self.setParam(nowdate, np.nanmean(lat), np.nanmean(sunz), np.nanmean(suna),
                      SatID, InstID, BandId, dem=dem)
        Coeffa, Coeffb, Coeffc  = self.corrCoeff()

        #大气校正
        corrdata = self.corrImage(radiance, Coeffa, Coeffb, Coeffc)

        print('第%s波段大气校正完成' %BandId)

        return corrdata

    def setParam(self, nowdate, lat, sunz, suna, SatID, InstID, BandId, dem=0.010):
        '''

        :param nowdate:
        :param BandId:
        :param radiance:
        :param metadata:
        :param dem:
        :return:
        '''

        self.set_geom(nowdate, sunz=sunz, suna=suna, satz=0, sata=0)
        self.set_atm(sLatitude=lat, nowdate=nowdate)
        self.set_aer()
        self.set_vis()
        self.set_altitude(dem=dem)

        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        SRFFile = os.path.join(EXEPATH, 'resp', 'FY', "FY.json")
        if not os.path.isfile(SRFFile) :
            raise Exception('风云卫星光谱响应文件不存在【%s】' %(SRFFile))
        SRF = json.load(open(SRFFile))

        minwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][0]
        maxwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][1]
        response = SRF[SatID][InstID]['B%d' %(BandId)]['SRF']
        self.set_resp(minwl, maxwl, response)


    def Calibration(self, filename, sdsname):

        # EV_250_Aggr.1KM_RefSB     1~4         # 0.47, 0.55, 0.65, 0.865,
        # EV_1KM_RefSB              5~19        # 1.38, 1.64, 2.13, 0.412,
        #                                         0.443, 0.49, 0.555, 0.67,
        #                                         0.709, 0.746, 0.865, 0.905,
        #                                         0.936, 0.94, 1.03,
        # EV_1KM_Emissive           20~23       # 3.796, 4.046, 7.233, 8.56,
        # EV_250_Aggr.1KM_Emissive  24~25       # 10.714, 11.948

        if 'EV_1KM_Emissive' in sdsname :
            data = readhdf(filename, sdsname)
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return self.calemiss(data, 10000.0/wave_length[np.arange(20, 24)-1])

        elif 'EV_250_Aggr.1KM_Emissive' in sdsname :
            data = readhdf(filename, sdsname)
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return self.calemiss(data, 10000.0/wave_length[np.arange(24, 26)-1])

        elif 'EV_250_Aggr.1KM_RefSB' in sdsname :
            data = readhdf(filename, sdsname)
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')

            return self.calref(data, coef[np.arange(1, 5)-1, :])

        elif 'EV_1KM_RefSB' in sdsname :
            data = readhdf(filename, sdsname)
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')

            return self.calref(data, coef[np.arange(5, 20)-1, :])

        else:
            tmpdata, sdsinfo = readhdf(filename, sdsname,dictsdsinfo={})
            data = tmpdata.copy()
            if 'Slope' in sdsinfo and 'Intercept' in sdsinfo:
                data = data * sdsinfo['Slope'] + sdsinfo['Intercept']

            if 'valid_range' in sdsinfo :
                data[(tmpdata < sdsinfo['valid_range'][0]) | (tmpdata > sdsinfo['valid_range'][1])] = np.nan

            return data

    def calref(self, dn, coef):
        ''' 可见光通道定标 '''
        ref = np.full_like(dn, fill_value=np.nan, dtype=np.float32)
        for i in range(dn.shape[0]) :
            ref[i] = coef[i, 0] + coef[i, 1] * dn[i] + coef[i, 2] * dn[i] * dn[i]

            ref[i, dn[i] == 65535] = np.nan
            ref[i, ref[i]<=0] = 0

        return ref

    def calemiss(self, dn, wavenum):
        ''' 红外波段定标 '''
        temp = dn * 0.01
        temp[(temp<=0) | (temp >= 600.0)] = np.nan

        bt = np.full_like(temp, fill_value=np.nan, dtype=np.float32)

        for i in range(temp.shape[0]) :
            bt[i] = self.planck_r2t(temp[i, :, :], wavenum[i])

        return bt

    def planck_r2t(self, rad, wn):
        '''
        普朗克函数：将辐射值转成亮温（K）

        Parameters
        ----------
        rad : numpy.narray
            mW/(m2.cm-1.sr)
        wn : float or numpy.narray
            wave number(cm^-1)

        Returns
        -------
            numpy.narray
            bright temperature
            units: K
        '''
        # 普朗克系数
        Radiation_C1 = 1.191042953E-5
        Radiation_C2 = 1.4387774

        bt = (Radiation_C2 * wn / np.log(Radiation_C1 * wn * wn * wn / (rad)+1.0))

        return bt

    def getresp(self, filename):
        # 对卫星中心提供的光谱响应函数解析

        data = np.loadtxt(filename, dtype=np.float32, skiprows=4)
        wavelength = data[:, 0] * 0.001 #
        response = data[:, 1]
        # wavelength = 10000.0 / wavenumber

        flag = response > 0.2
        wavelength = wavelength[flag]
        response = response[flag]

        minwl = np.ceil(np.nanmin(wavelength) / 0.0025) * 0.0025
        maxwl = np.floor(np.nanmax(wavelength) / 0.0025) * 0.0025

        fit = interpolate.interp1d(wavelength, response, kind= 'slinear')

        newwl = np.arange(minwl, maxwl+0.0025, 0.0025)
        minwl = np.nanmin(newwl)
        maxwl = np.nanmax(newwl)
        response = fit(newwl)

        return minwl, maxwl, response

