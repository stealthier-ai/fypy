# -*- coding:utf-8 -*-
'''
@Project     :

@File        :

@Modify Time :

@Author      : Lee

@Version     : 1.0

@Description :
绘制各种各样的色标卡
'''

import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont


def getColorList(cbfile, vmin=None, vmax=None):
    ''' 解析指定格式的colorbar文件 '''

    colorlist = {}
    cbset = {}
    # print(cbfile)
    if not os.path.isfile(cbfile) :
        raise Exception('文件不存在【%s】' %(cbfile))
        return None, cbset

    try:
        fp = open(cbfile, encoding='utf-8')
    except BaseException as e:
        fp = open(cbfile, encoding='gbk')

    lines = fp.readlines()

    line0 = lines[0].split(',')
    for key  in line0 :
        if not ':' in key :
            continue

        key = key.lower()
        keyvalue = key.split(':')[1].replace('\n', '')
        if 'cbtitle' in key :
            if len(keyvalue.replace(' ','')) == 0 :
                continue

            cbtitle = keyvalue
            if '(' in keyvalue :
                tls = keyvalue.split('(')
                cbtitle = tls[0]
                for item in tls[1:] :
                    cbtitle += '\n(' + item

            cbset['cbtitle'] = cbtitle
        elif 'cbdir' in key :
            cbdir = keyvalue.replace(' ', '')
            if len(cbdir) > 0 :
                cbset['cbdir'] = cbdir
        elif 'split' in key :
            split = keyvalue.replace(' ', '')
            if len(split) == 0 :
                split = False
            else:
                if split.lower() in ['true', 't'] :
                    split = True
            cbset['split'] = split
        elif 'row' in key :
            row = keyvalue.replace(' ', '')
            try:
                irow = int(row)
                if irow > 0 :
                    cbset['row'] = irow
            except BaseException :
                continue
        elif 'col' in key :
            col = keyvalue.replace(' ', '')
            try:
                icol = int(col)
                if icol > 0 :
                    cbset['col'] = icol
            except BaseException :
                continue
        elif 'vmin' in key :
            strvmin = keyvalue.replace(' ', '')
            try:
                vmin = float(strvmin)
                cbset['vmin'] = vmin
            except BaseException :
                continue
        elif 'vmax' in key :
            strvmax = keyvalue.replace(' ', '')
            try:
                vmax = float(strvmax)
                cbset['vmax'] = vmax
            except BaseException :
                continue

    # 渲染方式
    key = lines[1].split(':')[1].replace('\n', '')
    colorlist[key] = []
    for line in lines[2:] :
        line = line.replace('\n', '')
        item = line.split(',')
        if len(item) < 6 :
            continue
        # 6要素格式 ： value, r, g, b, alpha, label
        colorlist[key].append([float(item[0]),
                               (int(item[1]), int(item[2]), int(item[3]), int(item[4])),
                               item[5]])

    # 重构数据分段
    if vmin is not None and vmax is not None :
        for key in colorlist :
            if key in ['EXACT', 'DISCRETE'] :
                continue

            if vmin is None :
                vmin = colorlist[key][0][0]
            if vmax is None :
                vmax = colorlist[key][-1][0]

            cnt = 0
            for value, color, label in colorlist[key] :
                step = (vmax - vmin) / (len(colorlist[key]) -1)
                vvalue = vmin + step * cnt
                colorlist[key][cnt][0] = vvalue
                cnt  += 1

    return colorlist, cbset


class colorbar(object) :

    def __init__(self, cblist, width=1260, height=180,
                 cbtitle=None, titlesize=55,
                 ticks=None, ticklabels=None, tickfontsize=45,
                 vmin=None, vmax=None, orientation=None,
                 fmt='%.1f', font=None, FontTTF=None,
                 loc='center', split=False, cbset=None,
                 outname=None, **kwargs):
        '''
        loc : str
            tick的位置，'center' or 'centre', 'left', 'right'
                       'center' or 'centre', 'top', 'bottom'
        kwargs
        '''

        self.cblist = self.checkcblist(cblist)

        self.set_cbsize(width, height)

        if cbset is None : cbset={}

        self.set_vmin(vmin, cbset=cbset)
        self.set_vmax(vmax, cbset=cbset)

        # 设置图例方向 h:水平 v:垂直
        self.set_cbdir(orientation, cbset=cbset)

        # 设置图例块是否分隔
        self.set_cbsplit(split=split, cbset=cbset)

        # 设置标题以及字体、大小
        self.set_cbtitle(cbtitle=cbtitle, fontsize=titlesize, font=font, FontTTF=FontTTF)

        self.set_ticks(ticks=ticks, ticklabels=ticklabels, tickfontsize=tickfontsize,
                       loc=loc, fmt=fmt, font=font, FontTTF=FontTTF)

        self._colorkwargs = kwargs

        if outname is not None :
            self.save(outname)

    def save(self, outname, **kwargs):

        self.drawcb(**kwargs)

        # 保存colorbar为图片文件
        self.cbimg.save(outname)

    def show(self, **kwargs):

        self.drawcb(**kwargs)

        self.cbimg.show()

    def drawcb(self, alpha=0, **kwargs):

        # 创建画布
        self.cbimg = Image.new('RGBA', (int(self.cbwidth), int(self.cbheight)),
                               (255, 255, 255, alpha))
        self.cbdraw = ImageDraw.Draw(self.cbimg)

        cblist = self.cblist
        ticks = self.ticks
        ticklabels = self.ticklabels
        title = self.title

        if 'EXACT' in cblist:           # 分类渲染
            self.draw_exact_colorbar_img(cblist['EXACT'],
                                         ticks=ticks, title=title, **kwargs)
        elif 'DISCRETE' in cblist:      # 分级渲染
            self.interp = False
            if self.split :
                self.draw_exact_colorbar_img(cblist['DISCRETE'],
                                             ticks=ticks, title=title, **kwargs)
            else:
                self.draw_discrete_colorbar_img(cblist['DISCRETE'],
                                                ticks=ticks, ticklabels=ticklabels, title=title,
                                                **kwargs)
        elif 'INTERPOLATED' in cblist : # 连续渐变
            self.interp = True
            self.draw_interpolated_colorbar_img(cblist['INTERPOLATED'],
                                                ticks=ticks, ticklabels=ticklabels, title=title, **kwargs)

    def draw_exact_colorbar_img(self, cb, ticks=None, title=None, rows=None, cols=None,
                                outline=(0,0,0), width=2, **kwargs):
        '''
        绘制分类调色板，例如云分类、地表类型
        Parameters
        ----------
        cb
        ticks
        title
        outline
        width
        kwargs

        Returns
        -------

        '''

        if ticks is None :
            ticks = [label[2] for label in cb]

        colorbarcount = len(ticks)
        if colorbarcount == 0 :
            return

        if self.colorbardir :
            self.set_exact_cb_h(cb, colorbarcount, title, rows=rows, cols=cols, outline=outline, width=width)
        else:
            self.set_exact_cb_v(cb, colorbarcount, title, rows=rows, cols=cols, outline=outline, width=width)

    def draw_discrete_colorbar_img(self, cb, ticks=None, title=None, ticklabels=None, **kwargs):
        '''
        # 绘制分类调色板，例如海温、海风等分级显示
        :param cb:
        :param ticks:
        :param title:
        :param kwargs:
        :return:
        '''
        # 处理colorbar颜色模块

        if self.colorbardir : # 水平colorbar
            self.set_discrete_cb_h(cb, title, ticks, ticklabels)
        else:
            self.set_discrete_cb_v(cb, title, ticks, ticklabels)

    def draw_interpolated_colorbar_img(self, cb, ticks=None, ticklabels=None, title=None, **kwargs):

        # 处理colorbar颜色模块
        if self.colorbardir :
            self.set_interpolated_cb_h(cb, title, ticks, ticklabels)
        else:
            self.set_interpolated_cb_v(cb, title, ticks, ticklabels)

    def set_cbsize(self, width=1260, height=180):
        self.cbwidth = width
        self.cbheight = height

    def set_cbtitle(self, cbtitle, fontsize=None, font=None, FontTTF=None):
        self.title = cbtitle

        if font is None :
            self._titlefont = self.setfont(FontTTF, fontsize)
        else :
            self._titlefont = font

        if fontsize is not None :
            self.set_cbtitle_font(FontTTF=FontTTF, fontsize=fontsize)

    def set_cbtitle_font(self, FontTTF, fontsize):

        self.fontsize = fontsize
        self._titlefont = ImageFont.truetype(font=FontTTF, size=fontsize)

    def set_vmin(self, vmin=None, cbset=None):
        ''' 设置最小值 '''

        if vmin is None and 'vmin' in cbset :
            vmin = cbset['vmin']
        self.vmin = vmin

    def set_vmax(self, vmax=None, cbset=None):
        ''' 设置最大值 '''

        if vmax is None and 'vmax' in cbset :
            vmax = cbset['vmax']
        self.vmax = vmax

    def set_cbdir(self, cbdir, cbset=None):
        ''' 设置图例方向 '''

        if cbdir is None and 'cbdir' in cbset :
            cbdir = cbset['cbdir']

        if cbdir is None :
            self.colorbardir = True
        elif cbdir in ['h', 'horizontal'] :  # 水平
            self.colorbardir = True
        else:                                # 垂直
            self.colorbardir = False

    def set_cbsplit(self, split=None, cbset=None):

        if split is None and 'split' in cbset :
            split = cbset['split']

        self.split = split

    def set_ticks(self, ticks=None, ticklabels=None, tickfontsize=None,
                  loc=None, fmt=None, font=None, FontTTF=None):

        self.ticks = ticks
        if font is None :
            self._colorbarfont = self.setfont(FontTTF, tickfontsize)
        else :
            self._colorbarfont = font

        self._fmt = fmt
        self.loc = loc

        self.ticks = ticks
        self.ticklabels = ticklabels

    def set_ticklabels(self, ticklabels, fmt='%.1f'):

        self.ticklabels = ticklabels
        self._fmt = fmt

    def set_ticklabel_loc(self, loc):
        '''

        Parameters
        ----------
        loc : str
            tick的位置，'center' or 'centre', 'left', 'right'
                       'center' or 'centre', 'top', 'bottom'

        Returns
        -------

        '''

        self.loc = loc

    def set_ticklabel_font(self, FontTTF, fontsize):
        ''' 设置tick label的字体和字体大小 '''

        self.fontsize = fontsize

        self._colorbarfont = ImageFont.truetype(font=FontTTF, size=fontsize)

    def set_exact_cb_h(self, cb, colorbarcount, title, rows=1, cols=None, outline=(0,0,0), width=2):
        ''' 绘制分类横向colorbar '''

        # 设置Y向位置
        ratio_cb_offy = 0.1
        ratio_cb_y = 0.33

        # 设置X向位置
        ratio_space_x = 0.033
        ratio_cb_x = 0.93

        if isinstance(title, str) :
            title_left = 0
            title_top = 0
            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw,
                                                             title_left, title_top, title, drawflag=False)
        else:
            fontwidth = 0
            fontheight = 0

        # 如果未指定colorbar的行数，默认只一行
        if rows is None :
            rows = 1
            cols = colorbarcount

        #矩形间隔
        if colorbarcount <= 4 :
            # 矩形大小，不足4个，按照4个进行设置
            colorbarcount = 4
            ractangle_w = int((self.cbwidth-fontwidth) * 0.2)
            ractangle_h = int(self.cbheight * 0.5)
            top = int(self.cbheight * 0.15)
            left = int(((self.cbwidth-fontwidth) - ractangle_w*colorbarcount) / (colorbarcount+1))
            ranctangle_interval = int(((self.cbwidth-fontwidth) - ractangle_w*colorbarcount) / (colorbarcount+1))
        else:
            # 矩形大小
            ractangle_w = int((self.cbwidth-fontwidth) / colorbarcount * 0.7)
            ractangle_h = int(self.cbheight * 0.5)
            top = int(self.cbheight * 0.15)
            left = int(((self.cbwidth-fontwidth) - ractangle_w*colorbarcount) / (colorbarcount+1))
            ranctangle_interval = int(((self.cbwidth-fontwidth) - ractangle_w*colorbarcount) / (colorbarcount+1))

        for item in cb:
            label = None
            if len(item) == 3:
                val, rgba, label = item
            elif len(item) == 2 :
                val, rgba = item
            else:
                raise Exception('single colorbar error, '
                                'please check the cb["single"] format is '
                                '[value,[r, g, b], "label"] or [value,[r, g, b]]')

            #左上角的点
            x1, y1 = left, top
            #右下角的点
            x2 = x1 + ractangle_w
            y2 = y1 + ractangle_h

            left += ractangle_w + ranctangle_interval

            self.cbdraw.rectangle((x1, y1, x2, y2), fill=tuple(rgba), outline=outline, width=width)

            posx = int(x1)
            posy = int(y1 + ractangle_h)
            self.draw_ticks_h(self.cbdraw, posx, posy, val,
                              color_w=ractangle_w, color_h=ractangle_h, label=label)

        if isinstance(title, str) :
            title_left = int(left - ranctangle_interval*0.5)
            title_top = int(top + (ractangle_h - fontheight)*0.5)

            self._set_colorbar_title(self.cbdraw, title_left, title_top, title)

    def set_exact_cb_v(self, cb, colorbarcount, title, rows=None, cols=None, outline=(0,0,0), width=2):

        fontwidth = 0
        fontheight = 0
        top = int(self.cbheight * 0.02)
        left = int(self.cbwidth * 0.02)
        if isinstance(title, str) :
            title_left = int(left)
            title_top = int(top)

            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw, title_left, title_top, title, drawflag=True)


        # 矩形大小，不足4个，按照4个进行设置
        if colorbarcount <= 4 :
            colorbarcount = 4

        # 如果指定了行数，计算列数
        if rows is not None :
            cols = int(np.ceil(colorbarcount/rows))

        if cols is None :
            cols = 1

        # 根据指定的列数计算行数
        if cols is not None :
            rows = int(np.ceil(colorbarcount/cols))

        # 矩形大小
        ractangle_w = int(self.cbwidth / cols * 0.2)
        ractangle_h = int((self.cbheight - (top+fontheight)) / (rows * 4.0/3))

        # 两个图之间的间隔是图高度的33%
        ranctangle_interval = int(ractangle_h * 1.0 / 3)

        # 根据title和间隔进行调整colorbar的位置
        top = int((top+fontheight+ranctangle_interval/2))

        index = 0
        for j in range(cols) :
            for i in range(rows) :
                if index >= len(cb) :
                    break

                item = cb[index]
                index += 1

        # for item in cb:
                label = None
                if len(item) == 3:
                    val, rgba, label = item
                elif len(item) == 2 :
                    val, rgba = item
                else:
                    raise Exception('single colorbar error, '
                                    'please check the cb["single"] format is '
                                    '[value,[r, g, b], "label"] or [value,[r, g, b]]')

                #左上角的点
                x1 = left + j * int(self.cbwidth / cols)
                y1 = top  + i * (ractangle_h + ranctangle_interval)
                #右下角的点
                x2 = x1 + ractangle_w
                y2 = y1 + ractangle_h

                # top += ractangle_h + ranctangle_interval
                # 绘制矩形框
                self.cbdraw.rectangle((x1, y1, x2, y2), fill=tuple(rgba), outline=outline, width=width)

                # 计算ticklabel的起始xy位置
                posx = int((x1 + ractangle_w))
                posy = int(y1)

                # 绘制ticks
                self.draw_ticks_v(self.cbdraw, posx, posy, val,
                                  color_w=ractangle_w, color_h=ractangle_h, label=label)

    def set_discrete_cb_h(self, cb, title, ticks, ticklabels):
        # 设置Y向位置
        ratio_cb_offy = 0.1
        ratio_cb_y = 0.33

        # 设置X向位置
        ratio_space_x = 0.033
        ratio_cb_x = 0.93

        left = int(self.cbwidth * ratio_space_x)
        top = int(self.cbheight * ratio_cb_offy)

        if isinstance(title, str) :
            title_left = int(self.cbwidth * (ratio_space_x + ratio_cb_x))
            title_top = int(top)
            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw,
                                                             title_left, title_top, title, drawflag=False)
        else:
            fontwidth = 0
            fontheight = 0

        self._color_h = int(self.cbheight * ratio_cb_y)
        self._color_w = int((self.cbwidth-fontwidth) * ratio_cb_x)

        # 获取分级调色
        colorbar_data = self._mk_colorbar_data_h(cb, self._color_w, self._color_h)

        colorbar_data = self.get_rgb_color_data(colorbar_data, cb)

        colorbar_img = Image.fromarray(colorbar_data[:,:,:3], mode="RGB")
        # 将colorbar贴到图上
        self.cbimg.paste(colorbar_img, (left, top))

        # 设置ticks和ticklabels
        if ticks is None :
            ticks = []
            ticklabels = []
            for icount in range(0, len(cb)) :
                tick, rgba, label = cb[icount]
                ticks.append(tick)
                ticklabels.append(label)

        count= len(ticks)
        for icount in range(count) :
            tick = ticks[icount]
            label = ticklabels[icount]
            if len(label) == 0 :
                continue

            if tick > self.vmax or tick < self.vmin :
                continue

            posx = 1.0*icount / count * self._color_w + left
            posy = int(top + self._color_h*1.1)
            self.draw_ticks_h(self.cbdraw, posx, posy, tick,
                              color_w=self._color_w/count, color_h=self._color_h, label=label)

        if isinstance(title, str) :
            title_left = int(left*1.5 + self._color_w)
            title_top = int(top)

            self._set_colorbar_title(self.cbdraw, title_left, title_top, title)

    def set_discrete_cb_v(self, cb, title, ticks, ticklabels):
        # 垂直colorbar
        # 设置Y向位置
        ratio_cb_offy = 0.033  #  设置色带的上边空白百分比
        ratio_cb_y = 0.95   #  设置色带的高度百分比

        # 设置X向位置
        ratio_space_x = 0.02
        ratio_cb_x = 0.1  # 设置色带的宽度百分比

        top = int(self.cbheight * ratio_cb_offy)
        left = int(self.cbwidth * ratio_space_x)

        # 绘制title
        if isinstance(title, str) :
            title_left = int(left)
            title_top = int(self.cbheight * (ratio_cb_offy/2))
            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw, title_left, title_top, title)
            top = top + fontheight

        self._color_h = int((self.cbheight-top) * ratio_cb_y)
        self._color_w = int((self.cbwidth-left) * ratio_cb_x)

        # 绘制colorbar部分
        colorbar_data = self._mk_colorbar_data_v(cb, self._color_w, self._color_h)
        colorbar_data = self.get_rgb_color_data(colorbar_data, cb)
        colorbar_img = Image.fromarray(colorbar_data[:,:,:3], mode="RGB")
        self.cbimg.paste(colorbar_img, (left, top))

        # 设置ticks和ticklabels
        if ticks is None :
            ticks = []
            ticklabels = []
            for icount in range(0, len(cb)) :
                tick, rgba, label = cb[icount]
                ticks.append(tick)
                ticklabels.append(label)

        count= len(ticks)
        for icount in range(count) :
            tick = ticks[icount]
            label = ticklabels[icount]
            if len(label) == 0 :
                continue
            posx = int(left + self._color_w)
            posy = 1.0*icount / count * self._color_h + top
            self.draw_ticks_v(self.cbdraw, posx, posy, tick,
                              color_w=self._color_w, color_h=self._color_h/count, label=label)

    def set_interpolated_cb_h(self, cb, title, ticks, ticklabels):
        # 设置Y向位置
        ratio_cb_offy = 0.24
        ratio_cb_y = 0.26

        # 设置X向位置
        ratio_space_x = 0.033
        ratio_cb_x = 0.93

        top = int(self.cbheight * ratio_cb_offy)
        left = int(self.cbwidth * ratio_space_x)

        if isinstance(title, str) :
            title_left = int(self.cbwidth * (ratio_space_x + ratio_cb_x))
            title_top = int(top)
            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw, title_left, title_top, title, drawflag=False)
        else:
            fontwidth = 0
            fontheight = 0

        self._color_h = int(self.cbheight * ratio_cb_y)
        self._color_w = int((self.cbwidth-fontwidth) * ratio_cb_x)

        colorbar_data = self._mk_colorbar_data_h(cb, self._color_w, self._color_h)

        # 注意:由于设置的色带大小与数据范围重构后,导致色带大小变化,需要进行重置
        self._color_h, self._color_w = colorbar_data.shape

        # 绘制colorbar部分
        colorbar_data = self.get_rgb_color_data(colorbar_data, cb)
        colorbar_img = Image.fromarray(colorbar_data[:,:,:3], mode="RGB")
        self.cbimg.paste(colorbar_img, (left, top))

        ####################################################################
        # 设置刻度
        v_total = self.vmax - self.vmin
        # 设置ticks和ticklabels
        if ticks is None :
            ticks = []
            ticklabels = []
            for icount in range(0, len(cb)) :
                tick, rgba, label = cb[icount]
                if len(label) == 0 :
                    continue
                ticks.append(tick)
                ticklabels.append(label)

        if len(ticks) == 0 :
            ticks = [self.vmin + i * v_total*1.0 / 5 for i in range(6)]

        if len(ticklabels) == 0 or len(ticklabels) != len(ticks) :
            ticklabels = [self._fmt % (v, ) for v in ticks]

        count= len(ticks)
        for icount in range(count) :
            tick = ticks[icount]
            label = ticklabels[icount]
            if len(label) == 0 :
                continue

            if tick > self.vmax or tick < self.vmin :
                continue

            posx = 1.0 * (tick - self.vmin) / v_total * self._color_w + left
            posy = int(top + self._color_h*1.1)

            self.draw_ticks_h(self.cbdraw, posx, posy, tick,
                              color_w=self._color_w/v_total, color_h=self._color_h, label=label)

        if isinstance(title, str) :
            title_left = int(left*1.5 + self._color_w)
            title_top = int(top)
            self._set_colorbar_title(self.cbdraw, title_left, title_top, title)

    def set_interpolated_cb_v(self, cb, title, ticks, ticklabels):
        # 设置Y向位置
        ratio_cb_offy = 0.033  #  设置色带的上边空白百分比
        ratio_cb_y = 0.95   # 设置色带的高度百分比

        # 设置X向位置
        ratio_space_x = 0.02
        ratio_cb_x = 0.1  # 设置色带的宽度百分比

        top = int(self.cbheight * ratio_cb_offy)
        left = int(self.cbwidth * ratio_space_x)

        # 绘制title
        if isinstance(title, str) :
            title_left = int(left)
            title_top = int(self.cbheight * (ratio_cb_offy/2))
            fontwidth, fontheight = self._set_colorbar_title(self.cbdraw, title_left, title_top, title)
            top = top + fontheight

        self._color_h = int((self.cbheight-top) * ratio_cb_y)
        self._color_w = int((self.cbwidth-left) * ratio_cb_x)

        # 绘制colorbar部分
        colorbar_data = self._mk_colorbar_data_v(cb, self._color_w, self._color_h)
        # 注意:由于设置的色带大小与数据范围重构后,导致色带大小变化,需要进行重置
        self._color_h, self._color_w = colorbar_data.shape

        colorbar_data = self.get_rgb_color_data(colorbar_data, cb)
        colorbar_img = Image.fromarray(colorbar_data[:,:,:3], mode="RGB")
        self.cbimg.paste(colorbar_img, (left, top))

        v_total = self.vmax - self.vmin
        # 设置ticks和ticklabels
        if ticks is None :
            ticks = []
            ticklabels = []
            for icount in range(0, len(cb)) :
                tick, rgba, label = cb[icount]
                if len(label) == 0 :
                    continue
                ticks.append(tick)
                ticklabels.append(label)

        if len(ticks) == 0 :
            ticks = [self.vmin + i * v_total*1.0 / 5 for i in range(6)]

        if len(ticklabels) == 0 or len(ticklabels) != len(ticks) :
            ticklabels = [self._fmt % (v, ) for v in ticks]

        count= len(ticks)
        for icount in range(count) :
            tick = ticks[icount]
            label = ticklabels[icount]
            if len(label) == 0 :
                continue
            posx = int(left + self._color_w)
            posy = 1.0*(tick - self.vmin) / v_total * self._color_h + top
            self.draw_ticks_v(self.cbdraw, posx, posy, tick,
                              color_w=self._color_w, color_h=self._color_h/v_total, label=label)

    def _set_colorbar_title(self, draw, left, top, title=None, drawflag=True) :
        ''' 绘制colorbar的标题 '''

        if title is None :
            return None, None

        try:
            title = title.decode('utf8')
        except:
            pass

        if '\\n' in title :
            titles = title.split('\\n')
            fontwidth = 0
            fontheight = 0
            for key in titles :
                key = key.replace('\n', '')
                fontwidthtemp, fontheighttemp = self._titlefont.getsize(key)

                fontheight += fontheighttemp

                if fontwidthtemp > fontwidth :
                    fontwidth = fontwidthtemp
        else:
            fontwidth, fontheight = self._titlefont.getsize(title)

        if not drawflag :
            return fontwidth, fontheight

        if self.colorbardir :
            titles = title.split('\\n')
            count = 0
            for key in titles :
                key = key.replace('\n', '')
                fontwidthtemp, fontheighttemp = self._titlefont.getsize(key)
                posx = int(left+ fontwidth/2 - fontwidthtemp/2)
                posy = int(top+count*fontheighttemp*1.1)

                if posy < 0:
                    posy = 0

                self.cbwidth = posx

                if 'fontcolor' in self._colorkwargs :
                    draw.text((posx, posy), key, fill=tuple(self._colorkwargs['fontcolor']), font=self._titlefont)
                else:
                    draw.text((posx, posy), key, fill=(0, 0, 0), font=self._titlefont)

                count += 1
        else:
            titles = title.split('\\n')
            count = 0
            for key in titles :
                key = key.replace('\n', '')
                fontwidthtemp, fontheighttemp = self._titlefont.getsize(key)
                posx = int(left)
                posy = int(top+count*fontheighttemp)

                if posy < 0:
                    posy = 0

                if posx < 0:
                    posx = 0

                if 'fontcolor' in self._colorkwargs :
                    draw.text((posx, posy), key, fill=tuple(self._colorkwargs['fontcolor']), font=self._titlefont)
                else:
                    draw.text((posx, posy), key, fill=(0, 0, 0), font=self._titlefont)

                count += 1

        return fontwidth, fontheight

    def get_rgb_color_data(self, datas, cb):
        '''
        将数据值根据colorbar的设置，映射为r，g，b三色
        datas：二维数据值
        '''

        r = datas.copy()
        g = datas.copy()
        b = datas.copy()

        if self.interp :    # 连续渐变
            r_left_index = np.where(datas<=cb[0][0])
            g_left_index = r_left_index
            b_left_index = r_left_index
            r_right_index = np.where(datas>=cb[-1][0])
            g_right_index = r_right_index
            b_right_index = r_right_index

            for i in range(1, len(cb)):
                v_max = cb[i][0]
                v_min = cb[i-1][0]
                index = np.where((datas<v_max)&(datas>=v_min))

                #red ------
                c_max = cb[i][1][0]
                c_min = cb[i-1][1][0]
                kr = (c_max - c_min)*1.0 / (v_max - v_min)
                vr = kr*v_min - c_min
                r[index] *= kr
                r[index] -= vr

                #gren ------
                c_max = cb[i][1][1]
                c_min = cb[i-1][1][1]
                kg = (c_max - c_min)*1.0 / (v_max - v_min)
                vg = kg*v_min - c_min
                g[index] *= kg
                g[index] -= vg

                #blue ------
                c_max = cb[i][1][2]
                c_min = cb[i-1][1][2]
                kb = (c_max - c_min)*1.0 / (v_max - v_min)
                vb = kb*v_min - c_min
                b[index] *= kb
                b[index] -= vb

            last_index = np.where(datas==v_max)
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

            color_data = np.dstack([r, g, b])

            return color_data.astype(np.uint8)
        else:   # 分级
            r_left_index = np.where(datas<=cb[0][0])
            g_left_index = r_left_index
            b_left_index = r_left_index
            r_right_index = np.where(datas>=cb[-1][0])
            g_right_index = r_right_index
            b_right_index = r_right_index

            for i in range(1, len(cb)):
                v_min = cb[i-1][0]
                v_max = cb[i][0]
                index = (datas<=v_max)&(datas>v_min)

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

            return color_data.astype(np.uint8)

    # 水平方向从左到右创建colorbar数据
    def _mk_colorbar_data_h(self, cb, color_w, color_h):
        '''
        水平方向从左到右创建colorbar数据
        :param color_w:
        :param color_h:
        :return:
        '''
        pixel_items = []
        if self.interp : # 连续渐变
            self.vmin = cb[0][0]
            self.vmax = cb[-1][0]
            tval = self.vmax - self.vmin

            for i in range(1, len(cb)):
                # 值的间隔
                interval = cb[i][0] - cb[i-1][0]

                # 像素长度的百分比
                pixel_percent = (interval*1.0) / tval
                # 像素长度
                pixel_length = int(color_w * pixel_percent)
                if pixel_length == 0:
                    pixel_length = 1
                pixel_items.append((float(interval)/pixel_length,
                                    cb[i-1][0],
                                    cb[i][0],
                                    pixel_length))

            #每段数据列表
            count = 0
            broadcasting_list = []
            for item in pixel_items:
                count+=1
                # cnt = int((item[2] - item[1]) / item[0] )
                x = np.arange(item[1], item[2], item[0])
                x.shape = 1, len(x)
                #广播
                broadcasting_list.append(x.repeat(color_h, axis=0))
        else: # 分级
            # 值的间隔
            interval = cb[1][0] - cb[0][0]
            self.vmin = cb[0][0] - interval
            self.vmax = cb[-1][0]

            tval = self.vmax - self.vmin
            # 像素长度的百分比
            # pixel_percent = (interval*1.0) / tval
            pixel_percent = 1.0 / len(cb)
            # 像素长度
            pixel_length = int(color_w * pixel_percent)
            if pixel_length == 0:
                pixel_length = 1
            pixel_items.append((float(interval)/pixel_length,
                                self.vmin,
                                cb[0][0],
                                pixel_length))
            for i in range(1, len(cb)):
                # 值的间隔
                interval = cb[i][0] - cb[i-1][0]
                # 像素长度的百分比
                # pixel_percent = (interval*1.0) / tval
                pixel_percent = 1.0 / len(cb)
                # 像素长度
                pixel_length = int(color_w * pixel_percent)
                if pixel_length == 0:
                    pixel_length = 1
                pixel_items.append((float(interval)/pixel_length,
                                    cb[i-1][0],
                                    cb[i][0],
                                    pixel_length))

            #每段数据列表
            broadcasting_list = []
            for item in pixel_items:
                x = np.arange(item[1], item[2], item[0])
                x.shape = 1, len(x)
                #广播
                broadcasting_list.append(x.repeat(color_h, axis=0))

        return np.hstack(broadcasting_list)

    # 垂直方向从上到下创建colorbar数据
    def _mk_colorbar_data_v(self, cb, color_w, color_h):
        '''
        水平方向从左到右创建colorbar数据
        :param color_w:
        :param color_h:
        :return:
        '''
        pixel_items = []
        if self.interp : # 连续渐变
            for i in range(1, len(cb)):
                # 值的间隔
                interval = cb[i][0] - cb[i-1][0]
                self.vmin = cb[0][0]
                self.vmax = cb[-1][0]
                tval = self.vmax - self.vmin
                # 像素长度的百分比
                pixel_percent = (interval*1.0) / tval
                # 像素长度
                pixel_length = int(color_h * pixel_percent)
                if pixel_length == 0:
                    pixel_length = 1
                pixel_items.append((float(interval)/pixel_length,
                                    cb[i-1][0],
                                    cb[i][0],
                                    pixel_length))

            #每段数据列表
            broadcasting_list = []
            for item in pixel_items:
                x = np.arange(item[1], item[2], item[0])
                x.shape = len(x), 1
                #广播
                broadcasting_list.append(x.repeat(color_w, axis=1))
        else: # 分级
            # 值的间隔
            interval = cb[1][0] - cb[0][0]
            self.vmin = cb[0][0] - interval
            self.vmax = cb[-1][0]

            tval = self.vmax - self.vmin
            # 像素长度的百分比
            pixel_percent = (interval*1.0) / tval
            pixel_percent = 1.0 / len(cb)
            # 像素长度
            pixel_length = int(color_h * pixel_percent)
            if pixel_length == 0:
                pixel_length = 1
            pixel_items.append((float(interval)/pixel_length,
                                self.vmin,
                                cb[0][0],
                                pixel_length))
            for i in range(1, len(cb)):
                # 值的间隔
                interval = cb[i][0] - cb[i-1][0]
                # 像素长度的百分比
                pixel_percent = (interval*1.0) / tval
                pixel_percent = 1.0 / len(cb)
                # 像素长度
                pixel_length = int(color_h * pixel_percent)
                if pixel_length == 0:
                    pixel_length = 1
                pixel_items.append((float(interval)/pixel_length,
                                    cb[i-1][0],
                                    cb[i][0],
                                    pixel_length))

            #每段数据列表
            broadcasting_list = []
            for item in pixel_items:
                x = np.arange(item[1], item[2], item[0])
                x.shape = len(x), 1
                #广播
                broadcasting_list.append(x.repeat(color_w, axis=1))

        return np.vstack(broadcasting_list)

    def checkcblist(self, cblist):

        matchlist = {}
        for key in cblist :
            matchlist[key] = []
            for i in range(len(cblist[key])) :
                item = cblist[key][i]
                if len(item) == 3 and '$' in item[2] :
                    continue
                matchlist[key].append(item)

        return matchlist

    def setfont(self, FontTTF, fontsize):
        '''
        设置字体大小，如果没有传入字体大小，
        可根据图像的大小设置字体大小
        :param fontsize:
        :return:
        '''

        if fontsize is None :
            if self.colorbardir :
                if self.cbwidth < 300: fontsize = 5
                elif self.cbwidth < 700: fontsize = 9
                elif self.cbwidth < 1200: fontsize = 20
                elif self.cbwidth <= 2000: fontsize = 40
                elif self.cbwidth <= 4000: fontsize = 50
                elif self.cbwidth <= 6000: fontsize = 80
                elif self.cbwidth <= 13000: fontsize = 100
                elif  self.cbwidth <= 23000: fontsize = 110
                else: fontsize = 24
            else:
                if self.cbheight < 300: fontsize = 10
                elif self.cbheight < 600: fontsize = 18
                elif self.cbheight < 1200: fontsize = 40
                elif self.cbheight <= 2000: fontsize = 60
                elif self.cbheight <= 4000: fontsize = 80
                elif self.cbheight <= 6000: fontsize = 100
                elif self.cbheight <= 13000: fontsize = 110
                elif  self.cbheight <= 23000: fontsize = 120
                else: fontsize = 24
        self.fontsize = fontsize

        return ImageFont.truetype(font=FontTTF, size=fontsize)

    def draw_ticks_h(self, draw, posx, posy, val, color_w, color_h, label=None):
        if label is None :
            text = self._fmt % (val, )
        else:
            text = label
        try:
            text = text.decode('utf8')
        except:
            pass

        fontwidth, fontheight = self._colorbarfont.getsize(text)

        # 计算行列号
        if isinstance(self.loc, str) and self.loc in ['left'] :
            l = posx
        elif isinstance(self.loc, str) and self.loc in ['center', 'centre'] :
            l = posx + color_w*0.5
        else:
            l = posx + color_w

        if '&' in text :
            l = posx + color_w*0.5
            text = text.replace('&', '')

        titles = text.split('\\n')
        count = 0
        for key in titles :
            key = key.replace('\n', '')
            fontwidthtemp, fontheighttemp = self._titlefont.getsize(key)
            posx = int(l-fontwidthtemp/2)
            posy = int(posy+count*fontheighttemp)

            if posy < 0:
                posy = 0

            self.cbwidth = posx

            if 'fontcolor' in self._colorkwargs :
                draw.text((posx, posy), key, fill=tuple(self._colorkwargs['fontcolor']), font=self._colorbarfont)
            else:
                draw.text((posx, posy), key, fill=(0, 0, 0), font=self._colorbarfont)

            count += 1

        # posx = int(l - fontwidth/2)
        # posy = int(posy + color_h/8)
        # if posx < 0:
        #     posx = 0
        # if l + fontwidth/2 >= self.cbwidth :
        #     posx = self.cbwidth - fontwidth
        #
        # if 'fontcolor' in self._colorkwargs :
        #     draw.text((posx, posy), text, fill=tuple(self._colorkwargs['fontcolor']), font=self._colorbarfont)
        # else:
        #     draw.text((posx, posy), text, fill=(0,0,0), font=self._colorbarfont)

    def draw_ticks_v(self, draw, posx, posy, val, color_w, color_h, label=None):
        ''' 设置垂直colorbar的ticklabel'''
        if label is None :
            text = self._fmt % (val, )
        else:
            text = label
        try:
            text = text.decode('utf8')
        except:
            pass

        fontwidth, fontheight = self._colorbarfont.getsize(text)

        # 计算行列号
        if isinstance(self.loc, str) and self.loc in ['top'] :
            t = posy
        elif isinstance(self.loc, str) and self.loc in ['center', 'centre'] :
            t = posy + color_h*0.5
        else:
            t = posy + color_h

        if '&' in text :
            t = posy + color_h*0.5
            text = text.replace('&', '')

        posx = int(posx + color_w / 8)
        poxy = int(t - fontheight * 1/2)

        if posx + fontwidth >= self.cbwidth :
            posx = self.cbwidth - fontwidth

        if posx < 0:
            posx = 0

        if 'fontcolor' in self._colorkwargs :
            draw.text((posx, poxy), text, fill=tuple(self._colorkwargs['fontcolor']), font=self._colorbarfont)
        else:
            draw.text((posx, poxy), text, fill=(0,0,0), font=self._colorbarfont)



