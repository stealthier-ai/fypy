# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3config.py

@Modify Time :  2022/10/27 16:37   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from fypy import parm

exedir = os.path.abspath(list(parm.__path__)[0])

FontTTF = os.path.join(exedir, 'font', 'simsun.ttf')
if not os.path.isfile(FontTTF) :
    raise Exception('字体文件不存在【%s】' %(FontTTF))

FY4ProdInfo = {
    'FY4A':{
        'AGRI' : {
            'CIX': {
                'sdsname': ['Convective_Initiation'],
            },
            'CLM': {
                'sdsname': ['CLM'],
            },
            'CLP': {
                'sdsname': ['CLP'],
            },
            'CLT': {
                'sdsname': ['CLT'],
            },
            'CTH': {
                'sdsname': ['CTH'],
            },
            'CTP': {
                'sdsname': ['CTP'],
            },
            'CTT': {
                'sdsname': ['CTT'],
            },
            'DSD': {
                'sdsname': ['DST'],
            },
            'FHS': {
                'sdsname': ['FPA'],
            },
            'FOG': {
                'sdsname': ['FOG'],
            },
            'LST': {
                'sdsname': ['LST'],
            },
            'SST': {
                'sdsname': ['SST'],
            },
            'QPE': {
                'sdsname': ['Precipitation'],
            },
            'SNC': {
                'sdsname': ['SNC'],
            },
        },
        'LMI' : {
            'LMIE': {
                'sdsname': [''],
            },
            'LMIG': {
                'sdsname': [''],
            },
        },
        'GIIRS' : {
            'AVP': {
                'sdsname': ['AT_Prof'],
            },
        }
    }
}
