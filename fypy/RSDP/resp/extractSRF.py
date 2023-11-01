# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : extractSRF.py

@Modify Time :  2023/5/12 17:35   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import glob
import os
import sys
import numpy as np
import datetime

from scipy import interpolate
from fypy.tools import writejson

def getresp(wavelength, response):
    # 对提供的光谱响应函数解析

    # data = np.loadtxt(filename, dtype=np.float32, skiprows=4)
    # wavelength = data[:, 0] * 0.001 #
    # response = data[:, 1]
    # wavelength = 10000.0 / wavenumber

    flag = response > 0.1
    wavelength = np.array(wavelength[flag], dtype=np.float32)
    response = np.array(response[flag], dtype=np.float32)

    minwl = np.ceil(np.nanmin(wavelength) / 0.0025) * 0.0025
    maxwl = np.floor(np.nanmax(wavelength) / 0.0025) * 0.0025

    fit = interpolate.interp1d(wavelength, response, kind= 'slinear')

    newwl = np.arange(minwl, maxwl+0.0002, 0.0025)
    minwl = np.nanmin(newwl)
    maxwl = np.nanmax(newwl)
    response = fit(newwl)

    return minwl, maxwl, response


def extractFY():

    dict_info = {}
    filelist = glob.glob('./FY3D_MERSI/*.txt')
    for filename in filelist :
        data = np.loadtxt(filename, dtype=np.float32, skiprows=4)
        wavelength = data[:, 0] * 0.001 #
        response = data[:, 1]
        minwl, maxwl, response = getresp(wavelength, response)

        basename = os.path.basename(filename)
        namelist = basename.split('_')
        satid = namelist[0]
        instid = namelist[1]
        bandid = int(namelist[3][2:])

        if satid not in dict_info :
            dict_info[satid] = {}
        if instid not in dict_info[satid] :
            dict_info[satid][instid] = {}

        dict_info[satid][instid]['B%d' %(bandid)] = {}
        dict_info[satid][instid]['B%d' %(bandid)]['wl'] = [minwl, maxwl]
        dict_info[satid][instid]['B%d' %(bandid)]['SRF'] = list(response)

        print(minwl, maxwl)

    writejson('./FY/FY.json', dict_info)

def extractGF():
    import pandas as pd

    filename = r'D:\thinktanks\docs\高分卫星光谱响应\GF.xlsx'
    data = pd.read_excel(io=filename, sheet_name=None)

    dict_info = {}

    for key in data :
        keysplit = key.split('_')
        satid = keysplit[0]
        instid = keysplit[1]
        print(satid, instid)

        SRF = data[key].values

        wl = SRF[:, 0]
        resp = SRF[:,1:]

        cols, bands = resp.shape
        if satid not in dict_info :
            dict_info[satid] = {}
        if instid not in dict_info[satid] :
            dict_info[satid][instid] = {}

        for bandid in range(bands) :
            minwl, maxwl, response = getresp(wl*0.001, resp[:,bandid])

            dict_info[satid][instid]['B%d' %(bandid+1)] = {}
            dict_info[satid][instid]['B%d' %(bandid+1)]['SRF'] = list(response)
            dict_info[satid][instid]['B%d' %(bandid+1)]['wl'] = [minwl, maxwl]
            print(minwl, maxwl)

    writejson('./resp/GF/GF.json', dict_info)


def extractSentinel():

    dict_info = {}

    data = PredefinedWavelengths.__dict__

    for key in data :
        if isinstance(data[key], tuple) and len(data[key]) == 4:
            keysplit = key.split('_')
            if len(keysplit) == 3 :
                satid = keysplit[0]
                instid = keysplit[1]
                bandid = keysplit[2]
            else:
                satid = keysplit[2]
                instid = keysplit[1]
                bandid = keysplit[3]

            SRF = data[key]
            wv_id = SRF[0]
            if wv_id <= 0:
                continue

            # if satid not in ['S2A', 'S2B', 'S3A', 'S3B'] :
            #     continue

            if satid not in ['LANDSAT'] :
                continue

            wl = SRF[1:3]
            response = SRF[3]

            if satid not in dict_info :
                dict_info[satid] = {}
            if instid not in dict_info[satid] :
                dict_info[satid][instid] = {}

            try:
                bandname = 'B%d' %(int(bandid))
            except BaseException :
                bandname = bandid

            dict_info[satid][instid][bandname] = {}
            dict_info[satid][instid][bandname]['SRF'] = list(response)
            dict_info[satid][instid][bandname]['wl'] = wl

    writejson('./Landsat/Landsat.json', dict_info)

from Py6S import *
if __name__ == '__main__':

    # extractSentinel()

    extractFY()



