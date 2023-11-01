# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3L1Pro.py

@Modify Time :  2023/8/11 17:05   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from fypy.tools import readhdf, readhdf_fileinfo


class fy3L1Pro(object) :
    ''' 对FY3 MERSI L1数据进行辐射定标 '''

    def __init__(self, satID, instID):
        self.satID = satID
        self.instID = instID

    def calibration(self, filename, sdsname):
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
            return self.calref(data, coef[np.arange(1, 5)-1, :])
        elif 'EV_1KM_RefSB' in sdsname :
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')
            return self.calref(data, coef[np.arange(5, 20)-1, :])
        elif 'EV_1KM_Emissive' in sdsname :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return self.calemiss(data, 10000.0/wave_length[np.arange(20, 24)-1])
        elif 'EV_250_Aggr.1KM_Emissive' in sdsname :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            return self.calemiss(data, 10000.0/wave_length[np.arange(24, 26)-1])
        elif sdsname in ['EV_250_RefSB_b1', 'EV_250_RefSB_b2',
                         'EV_250_RefSB_b3', 'EV_250_RefSB_b4'] :
            coef = readhdf(filename, '/Calibration/VIS_Cal_Coeff')
            bandid = int(sdsname[-1])
            return self.calref(data, [coef[bandid-1, :]])
        elif sdsname in ['EV_250_Emissive_b24', 'EV_250_Emissive_b25'] :
            fileinfo = readhdf_fileinfo(filename)
            wave_length = fileinfo['Effect_Center_WaveLength']
            bandid = int(sdsname[-2:])
            return self.calemiss(data, [10000.0/wave_length[bandid-1]])
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

