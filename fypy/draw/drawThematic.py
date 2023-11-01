# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : drawThematic.py

@Modify Time :  2022/12/8 17:03   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from PIL import Image, ImageDraw, ImageFont
from fypy.draw.colorbar import getColorList

from fypy import parm

exedir = os.path.abspath(list(parm.__path__)[0])

FontTTF = os.path.join(exedir, 'font', 'simsun.ttf')
if not os.path.isfile(FontTTF) :
    raise Exception('字体文件不存在【%s】' %(FontTTF))


class drawThematic() :

    def __init__(self, data, cbfile=None, prodid=None, vmin=None, vmax=None):

        if cbfile is None :
            if prodid is None :
                raise Exception('请提供渲染方案或者产品ID')
            cbfile = os.path.join(exedir, 'cb', 'colorbar_%s.txt' %(prodid))

        srcdata = np.array(data)
        if len(srcdata.shape) == 2 :
            self.HEIGHT, self.WIDTH = srcdata.shape
            self.LEVEL = 1
        elif len(srcdata.shape) == 3 :
            self.LEVEL, self.HEIGHT, self.WIDTH = srcdata.shape
        else:
            raise Exception('当前仅支持2D或3D的数据')

        # 获取渲染方案
        colorlist, cbset = getColorList(cbfile, vmin=vmin, vmax=vmax)

        rgbarray = self.to_rgb(srcdata, colorlist)

        if rgbarray.shape[2] == 3:
            self.img = Image.fromarray(rgbarray, mode="RGB")
        elif rgbarray.shape[2] == 4:
            self.img = Image.fromarray(rgbarray, mode="RGBA")

    def to_rgb(self, srcdata, colorlist):
        '''
        将数据值根据colorbar的设置，映射为r，g，b三色

        Parameters
        ----------
        srcdata : numpy.array
            二维数据值
        colorlist : list
            渲染方案

        Returns
        -------
            [r, g, b]
        '''

        r = srcdata.copy()
        g = srcdata.copy()
        b = srcdata.copy()

        if 'INTERPOLATED' in colorlist :
            interp = True
            cb = colorlist['INTERPOLATED']
        elif 'DISCRETE' in colorlist :
            interp = False
            cb = colorlist['DISCRETE']
        elif 'EXACT' in colorlist :
            interp = False
            cb = colorlist['EXACT']
        else:
            raise Exception('当前仅支持连续渐变“INTERPOLATED“、分级渲染”DISCRETE“、分类渲染”EXACT”')

        if interp :    # 连续渐变
            r_left_index = np.where(srcdata <= cb[0][0])
            g_left_index = r_left_index
            b_left_index = r_left_index
            r_right_index = np.where(srcdata >= cb[-1][0])
            g_right_index = r_right_index
            b_right_index = r_right_index

            for i in range(1, len(cb)):
                v_max = cb[i][0]
                v_min = cb[i-1][0]
                index = np.where((srcdata < v_max) & (srcdata >= v_min))

                # 红
                c_max = cb[i][1][0]
                c_min = cb[i-1][1][0]
                kr = (c_max - c_min)*1.0 / (v_max - v_min)
                vr = kr*v_min - c_min
                r[index] *= kr
                r[index] -= vr

                # 绿
                c_max = cb[i][1][1]
                c_min = cb[i-1][1][1]
                kg = (c_max - c_min)*1.0 / (v_max - v_min)
                vg = kg*v_min - c_min
                g[index] *= kg
                g[index] -= vg

                # 蓝
                c_max = cb[i][1][2]
                c_min = cb[i-1][1][2]
                kb = (c_max - c_min)*1.0 / (v_max - v_min)
                vb = kb*v_min - c_min
                b[index] *= kb
                b[index] -= vb

            last_index = np.where(srcdata == v_max)
            r[last_index] *= kr
            r[last_index] -= vr
            g[last_index] *= kg
            g[last_index] -= vg
            b[last_index] *= kb
            b[last_index] -= vb

            #大于最大值的使用最大值的颜色，小于最小值的数据使用最小值的颜色
            r[r_left_index]  = cb[0][1][0]
            g[g_left_index]  = cb[0][1][1]
            b[b_left_index]  = cb[0][1][2]
            r[r_right_index] = cb[-1][1][0]
            g[g_right_index] = cb[-1][1][1]
            b[b_right_index] = cb[-1][1][2]

            # alpha = np.full_like(r, fill_value=cb[0][1][3])
            color_data = np.dstack([r, g, b])

            return color_data.astype(np.int8)
        else:   # 分级
            r_left_index = np.where(srcdata <= cb[0][0])
            g_left_index = r_left_index
            b_left_index = r_left_index
            r_right_index = np.where(srcdata >= cb[-1][0])
            g_right_index = r_right_index
            b_right_index = r_right_index

            for i in range(1, len(cb)):
                v_min = cb[i-1][0]
                v_max = cb[i][0]
                index = (srcdata <= v_max) & (srcdata > v_min)

                r[index] = cb[i][1][0]
                g[index] = cb[i][1][1]
                b[index] = cb[i][1][2]

            #大于最大值的使用最大值的颜色，小于最小值的数据使用最小值的颜色
            r[r_left_index]  = cb[0][1][0]
            g[g_left_index]  = cb[0][1][1]
            b[b_left_index]  = cb[0][1][2]
            r[r_right_index] = cb[-1][1][0]
            g[g_right_index] = cb[-1][1][1]
            b[b_right_index] = cb[-1][1][2]

            #rgb三色拼合
            color_data = np.dstack([r, g, b])

            return color_data.astype(np.int8)

    def text(self, x:float, y:float, s : str, fontsize=20, fontcolor=(255, 250, 110)):
        '''
        在图上绘制文本

        Parameters
        ----------
        x : float
            在X方向位置偏移百分比（0.0 - 1.0）
        y : float
            在Y方向位置偏移百分比（0.0 - 1.0）
        s ：string
            所需绘制的文本
        fontcolor： tuple or string (optional)
            文本的颜色，可设置RGB或者色标的标准名，.eg. "red"
        fontsize : int
            字体大小

        Returns
        -------

        '''
        x = np.int(self.WIDTH * x)
        y = np.int(self.HEIGHT * y)

        font = ImageFont.truetype(FontTTF, fontsize)

        draw = ImageDraw.Draw(self.img)
        draw.text((x, y), s, fill=fontcolor, font=font)

    def save(self, outname):
        ''' 保存为图片 '''

        self.img.save(outname)

