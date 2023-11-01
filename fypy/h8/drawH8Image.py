# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : drawH8Image.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
from PIL import Image


def drawH8TrueColor(outname, vis650, vis550, vis450):
    '''
    绘制葵花数据真彩图
    
    Parameters
    ----------
    outname : string
        输出文件名，PNG或者JPG
    vis650 : numpy.array-2D
        650um通道的反射率，范围在0-1
    vis550 : numpy.array-2D
        550um通道的反射率，范围在0-1
    vis450 : numpy.array-2D
        450um通道的反射率，范围在0-1

    Returns
    -------

    '''

    r = vis650 * 255
    g = vis550 * 255
    b = vis450 * 255

    r[r<0] = 0
    g[g<0] = 0
    b[b<0] = 0

    r[r>255] = 255
    g[g>255] = 255
    b[b>255] = 255

    rgbArray = np.zeros((r.shape[0], r.shape[1], 4), dtype=np.float64)
    rgbArray[..., 0] = r
    rgbArray[..., 1] = g
    rgbArray[..., 2] = b
    rgbArray[..., 3] = 255

    img = Image.fromarray(rgbArray.astype(np.uint8))
    img.save(outname, quality=90)


    # im = Image.merge('RGB', (Image.fromarray(np.array(r, dtype=np.uint8), "L"),
    #                          Image.fromarray(np.array(g, dtype=np.uint8), "L"),
    #                          Image.fromarray(np.array(b, dtype=np.uint8), "L")))
    # im.save(outname, quality=95)
