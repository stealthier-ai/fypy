# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4searchtable.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee

@Version     : 1.0

@Description :

'''
import os

import sys
import numpy as np
import datetime


TWOPI = 6.28318530717958648
DPAI = 6.28318530717958648

deg2rad = np.pi / 180.0
rad2deg = 180.0 / np.pi
FILLVALUE = -999


# 北边界纬度/º	80.56672132
# 南边界纬度/º	-80.56672132
# 东边界经度/º	-174.71662309
# 西边界经度/º	24.11662309


# 分辨率(km)   行数    列数
#   16        687    687
#   8        1374   1374
#   4        2748   2748
#   2        5496   5496
#   1       10992  10992
# 0.5       21984  21984
# 0.25      43968  43968

class fy4searchtable() :

    def __init__(self, subpoint, resolution):

        self.ea = 6378.137          # 地球长半轴
        self.eb = 6356.7523         # 地球短半轴
        self.H = 42164.0        # 地心到卫星质心的距离

        self.subpoint = subpoint
        self.resolution = resolution

        if resolution == 0.0025 :  # 250米
            self.coff = 21983.5
            self.cfac = 163730199
            self.loff = 21983.5
            self.lfac = 163730199
            self.rowmax = 43968
            self.colmax = 43968
        elif resolution == 0.005 :
            self.coff = 10991.5
            self.cfac = 81865099
            self.loff = 10991.5
            self.lfac = 81865099
            self.rowmax = 21984
            self.colmax = 21984
        elif resolution == 0.01 :
            self.coff = 5495.5
            self.cfac = 40932549
            self.loff = 5495.5
            self.lfac = 40932549
            self.rowmax = 10992
            self.colmax = 10992
        elif resolution == 0.02 :
            self.coff = 2747.5
            self.cfac = 20466274
            self.loff = 2747.5
            self.lfac = 20466274
            self.rowmax = 5496
            self.colmax = 5496
        elif resolution == 0.04 :
            self.coff = 1373.5      # 1375.5
            self.cfac = 10233137    # 10233137
            self.loff = 1373.5      # 1375.5
            self.lfac = 10233137    # 10233137
            self.rowmax = 2748
            self.colmax = 2748
        else:
            raise Exception("resolution error! please input [0.005, 0.01, 0.02, 0.04]")


        self.FLAT = 0.00335281317789691
        self.AE = 6378.137
        self.E2 = 0.0066943800699785


    def xytolatlon(self, x, y):
        '''
        行列号转经纬度

        Parameters
        ----------
        x:  numpy.array
            列号
        y:  numpy.array
            行号

        Returns
        -------
            numpy.array
            返回纬度、经度
        '''

        row = np.array(y)
        col = np.array(x)

        fillflag = (row < 0) | (row >= self.rowmax) | (col < 0) | (col >= self.colmax)

        x = deg2rad * (col - self.coff) / (2**-16 * self.cfac)
        y = deg2rad * (row - self.loff) / (2**-16 * self.lfac)

        sd = (self.H * np.cos(x) * np.cos(y)) * (self.H * np.cos(x) * np.cos(y)) \
             - (np.cos(y) * np.cos(y) + (self.ea * self.ea) / (self.eb * self.eb) \
                * np.sin(y) * np.sin(y)) * ((self.H * self.H) - self.ea * self.ea)

        flag = sd < 0

        sd[sd>=0] = np.sqrt(sd[sd>=0])

        sn = (self.H * np.cos(x) * np.cos(y) - sd) / \
             (np.cos(y) * np.cos(y) + (self.ea * self.ea) / (self.eb * self.eb) * np.sin(y) * np.sin(y))

        S1 = self.H - (sn * np.cos(x) * np.cos(y))
        S2 = sn * np.sin(x) * np.cos(y)
        S3 = -sn * np.sin(y)
        Sxy = np.sqrt(S1 * S1 + S2 * S2)

        lon = rad2deg * np.arctan2(S2 , S1) + self.subpoint
        lat = rad2deg * np.arctan((self.ea * self.ea) / (self.eb * self.eb) * S3 / Sxy)

        lon[lon > 180] -= 360.0
        lon[lon < -180] += 360.0

        lon[flag] = FILLVALUE
        lat[flag] = FILLVALUE

        lon[fillflag] = FILLVALUE
        lat[fillflag] = FILLVALUE

        lat = np.array(lat, dtype=np.float32)
        lon = np.array(lon, dtype=np.float32)

        return  lat, lon


    def latlon2xy(self, lat, lon):
        '''
        经纬度转行列号

        Parameters
        ----------
        lat : numpy.array
            纬度， -90.0 ~ 90.0
        lon : numpy.array
            经度， -180.0 ~ 180.0

        Returns
        -------
             numpy.array
             返回行号、列号
        '''

        lat1 = lat.copy()
        lon1 = lon.copy()

        lon1[lon1 > 180] -= 360.0

        lon1 = lon * deg2rad
        lat1 = lat * deg2rad

        phi_e = np.arctan(self.eb ** 2 * np.tan(lat1) / self.ea ** 2)
        re = self.eb / np.sqrt(1-(self.ea**2 - self.eb**2) * np.cos(phi_e)**2/self.ea**2)
        r1 = self.H - re * np.cos(phi_e) * np.cos(lon1 - self.subpoint * deg2rad)
        r2 = -re * np.cos(phi_e) * np.sin(lon1 - self.subpoint * deg2rad)
        r3 = re * np.sin(phi_e)
        rn = np.sqrt(r1**2 + r2**2 + r3**2)
        x = np.arctan2(-r2, r1)*180/np.pi
        # x = np.arctan(-r2/r1)*180/np.pi
        y = np.arcsin(-r3/rn)*180/np.pi
        col = self.coff + x * 2**(-16) * self.cfac
        line = self.loff + y * 2**(-16) * self.lfac

        # 对外太空进行判断
        tmp = r1 * (r1 - self.H) + r2 **2 + r3**2
        flag = tmp > 0

        line = np.array(line+0.5, dtype=np.int32)
        col = np.array(col+0.5, dtype=np.int32)

        line[(line<0) | (line >= self.rowmax)] = FILLVALUE
        col[(col<0) | (col>=self.colmax)]  = FILLVALUE

        return col, line


    def CalSatZenith(self, lat, lon) :
        lat1 = np.array(lat, dtype=np.float32)
        lon1 = np.array(lon, dtype=np.float32)
        fillflag = (lat1 < -90) | (lat1 > 90) | (lon1 < -180) | (lon1 > 180)

        lat1 = lat1 * deg2rad
        lon1 = lon1 * deg2rad
        point = self.subpoint * deg2rad

        c_lat = np.arctan(0.993243 * np.tan(lat1))
        rl = 6356.7523 / np.sqrt(1 - 0.00675701 * np.cos(c_lat) * np.cos(c_lat))
        r1 = self.H - rl * np.cos(c_lat) * np.cos(lon1 - point)
        r2 = -rl * np.cos(c_lat) * np.sin(lon1 - point)

        r3 = rl * np.sin(c_lat)
        rn = np.sqrt(r1 * r1 + r2 * r2 + r3 * r3)

        r4 = np.arccos((rn * rn + self.H * self.H - rl * rl) / (2 * rn * self.H))

        r5 = np.arcsin(np.sin(r4) * self.H / 6356.7523) * rad2deg

        r5[fillflag] = FILLVALUE

        return r5


    def CalSatAzimuth(self, lat, lon):
        '''
        计算卫星方位角
        :param lat1: array_like,
                    unit: degree
        :param lon1: array_like,
                    unit: degree

        :return: satellite azimuth , array_like:
                    unit:degree
        '''
        lat1 = np.array(lat, dtype=np.float32)
        lon1 = np.array(lon, dtype=np.float32)
        fillflag = (lat1 < -90) | (lat1 > 90) | (lon1 < -180) | (lon1 > 180)

        lat1 = lat1 * deg2rad
        lon1 = lon1 * deg2rad
        point = self.subpoint * deg2rad

        heigh = self.H * 1.0E3      # km -> m

        RX = self._Zocgef(lat1, lon1, heigh)

        RSat = self._Zocgef(0.0, point, heigh)

        Azi = self._ICGSSAT(lat1, lon1, RSat, RX)

        sata = Azi * rad2deg

        flag = sata >= 180
        sata[flag] -= 180.0
        sata[~flag] += 180.0

        sata[fillflag] = FILLVALUE

        return sata


    def CalSunZenith(self, lat, lon, nowtime):
        '''
        计算太阳天顶角
        :param lat: array_like,
                    unit: degree
        :param lon: array_like,
                    unit: degree
        :param nowtime: datetime
        :return: sun zenith , array_like:
                    unit:degree
        '''
        lat1 = np.array(lat, dtype=np.float32)
        lon1 = np.array(lon, dtype=np.float32)
        fillflag = (lat1 < -90) | (lat1 > 90) | (lon1 < -180) | (lon1 > 180)

        lat1 = lat1 * deg2rad
        lon1 = lon1 * deg2rad

        sttime = (nowtime - datetime.datetime.strptime('%s-01-01' % (nowtime.strftime('%Y')), '%Y-%m-%d')).total_seconds()

        PHIG,  DPHIDT, CHISUN,  EPSLON = self._EARTH(sttime)

        RSOL = 149600881.0
        COSX = np.cos(CHISUN)
        SINX = np.sin(CHISUN)
        COSE = np.cos(EPSLON)
        SINE = np.sin(EPSLON)

        RX = [RSOL*COSX,  RSOL*(SINX*COSE),  RSOL*(SINX*SINE)]

        TH = np.fmod(lon1 + PHIG, TWOPI)
        cTH = np.cos(TH)
        sTH = np.sin(TH)

        Azgs, Elgs = self._ICGS(RX, sTH, cTH, lat1)

        # SunAzi = Azgs * rad2deg
        SunZen = ((np.pi / 2.0) - Elgs) * rad2deg

        SunZen[fillflag] = FILLVALUE

        return SunZen


    def CalSunAzimuth(self, lat, lon, nowtime):

        lat1 = np.array(lat, dtype=np.float32)
        lon1 = np.array(lon, dtype=np.float32)
        fillflag = (lat1 < -90) | (lat1 > 90) | (lon1 < -180) | (lon1 > 180)

        lat1 = lat1 * deg2rad
        lon1 = lon1 * deg2rad

        sttime = (nowtime - datetime.datetime.strptime('%s-01-01' % (nowtime.strftime('%Y')), '%Y-%m-%d')).total_seconds()

        PHIG,  DPHIDT, CHISUN,  EPSLON = self._EARTH(sttime)

        RSOL = 149600881.0
        COSX = np.cos(CHISUN)
        SINX = np.sin(CHISUN)
        COSE = np.cos(EPSLON)
        SINE = np.sin(EPSLON)

        RX = [RSOL*COSX,  RSOL*(SINX*COSE),  RSOL*(SINX*SINE)]

        TH = np.fmod(lon1 + PHIG, TWOPI)
        cTH = np.cos(TH)
        sTH = np.sin(TH)

        Azgs, Elgs = self._ICGS(RX, sTH, cTH, lat1)

        SunAzi = Azgs * rad2deg
        # SunZen = ((np.pi / 2.0) - Elgs) * rad2deg

        SunAzi[fillflag] = FILLVALUE

        return SunAzi


    def CalSunGL(self, satz, sunz, rela):
        '''

        :param satz:
        :param sunz:
        :param rela:
        :return:
        '''
        fillflag = (satz < -360) | (satz > 360) | (sunz < -360) | (sunz > 360) | (rela < -360) | (rela > 360)

        data = np.cos(satz * deg2rad) * np.cos(sunz * deg2rad) \
               + np.sin(satz * deg2rad) * np.sin(sunz * deg2rad) * np.cos(rela * deg2rad)
        # cos_solzen * cos_satzen + sin_solzen * sin_satzen * cos_relaz

        data[data > 1.0] = 0.0
        SunGL = np.arccos(data) * rad2deg
        SunGL[fillflag] = -999.0

        return SunGL


    def CalRelativeAzimuth(self, sata, suna):
        fillflag = (sata > 360) | (sata < 0) | (suna > 360) | (suna < 0)

        RelativeAzi = np.fabs(sata - suna) + 180.0 # ???? + 180

        RelativeAzi[fillflag] = -999.0

        return RelativeAzi

    def _ICGSSAT(self, lat, lon, RSat, RX):

        HGT = 50.0 / 6378137.0
        Sla = np.sin(lat)
        Cla = np.cos(lat)
        Sth = np.sin(lon)
        Cth = np.cos(lon)

        Rhx = RSat[0] - RX[0]
        Rhy = RSat[1] - RX[1]
        Rhz = RSat[2] - RX[2]
        Rht = Cth * Rhx + Sth * Rhy
        Rxgs = Cth * Rhy - Sth * Rhx
        Rygs = Cla * Rhz - Sla * Rht
        Rzgs = Cla * Rht + Sla * Rhz

        Azgs = np.arctan2(Rxgs, Rygs)
        Azgs[Azgs < 0.0] = Azgs[Azgs < 0.0] + TWOPI
        Azgs[Azgs < np.pi] = Azgs[Azgs < np.pi] + np.pi
        Azgs[Azgs >= np.pi] = Azgs[Azgs >= np.pi] - np.pi

        Elgs = np.arctan(Rzgs / np.sqrt(Rxgs * Rxgs + Rygs * Rygs))

        return Elgs


    def _Zocgef(self, lat, lon, height):

        Xge = [] # np.zeros(3, dtype=np.float64)
        E = self.FLAT * (2.0 - self.FLAT)
        ENP = (self.AE * 1000) / np.sqrt(1.0 - E * np.sin(lat) * E * np.sin(lat))
        Xge.append( (ENP + height) * np.cos(lat) * np.cos(lon) / 1000.)

        Xge.append((ENP + height) * np.cos(lat) * np.sin(lon) / 1000.)

        Xge.append((ENP * (1.0 - E) + height) * np.sin(lat) / 1000.)

        return Xge


    def _ICGS(self, RX, Sth, Cth, lat):
        HGT = 50.0 / 6378137.0
        Sla = np.sin(lat)
        Cla = np.cos(lat)
        RE2 = self.AE / np.sqrt(1.0 - self.E2 * Sla * Sla)
        Gxic = (RE2 + HGT) * Cla * Cth

        Gyic = (RE2 + HGT) * Cla * Sth
        Gzic = (RE2 * (1.0 - self.E2) + HGT) * Sla
        Rhx = RX[0] - Gxic
        Rhy = RX[1] - Gyic
        Rhz = RX[2] - Gzic
        Rht = Cth * Rhx + Sth * Rhy
        Rxgs = Cth * Rhy - Sth * Rhx
        Rygs = Cla * Rhz - Sla * Rht
        Rzgs = Cla * Rht + Sla * Rhz

        Azgs = np.arctan2(Rxgs, Rygs)

        Azgs[Azgs < 0.0] = Azgs[Azgs < 0.0] + TWOPI
        Elgs = np.arctan(Rzgs / np.sqrt(Rxgs * Rxgs + Rygs * Rygs))

        return Azgs, Elgs

    def _EARTH(self, stime):
        PHIGRN = [ 1.739935890, 628.33195070, .000006754947810]
        CHIMS = [4.881633570, 628.3319487050, 1.4468088e-10]

        MEARTH = [6.256581470, 628.3019468590, -7.1676766e-11]
        OMEGA = [4.523636080, -33.7571619360, 9.9285594E-10]
        EPSLON = [.40936890950, - 2.27232172336e-4]
        KAPPA = .0335020
        DELPSI = [8.35465e-5,  .61712e-5]
        DELEPS = [4.46513e-5, .26771e-5]
        ZERO80 = 29219.50
        CENTCT = 36525.0

        DDATE = ZERO80 + stime / 86400.0
        IDATE2 = 2. * DDATE

        if IDATE2 % 2 == 0 :
            IDATE2 -= 1
        # IDATE2[(IDATE2 % 2) == 0] = IDATE2[(IDATE2 % 2) == 0] - 1
        DATE = IDATE2 * .50

        if DDATE - DATE > .50 :
            DATE += 1
        # DATE[DDATE - DATE > .50] = DATE[DDATE - DATE > .50] + 1.0
        EPTIME = (DATE - ZERO80) * 86400.0
        CENTNO = DATE / CENTCT

        ALPHA0 = self._ZSODSRadto2Pai(PHIGRN[0] + self._ZSODSRadto2Pai(CENTNO * (PHIGRN[1] + CENTNO * PHIGRN[2])))

        DADT = (TWOPI + (PHIGRN[1] + 2. * CENTNO * PHIGRN[2]) / CENTCT) / 86400.

        CHIMO = self._ZSODSRadto2Pai(CHIMS[0] + self._ZSODSRadto2Pai(CENTNO * (CHIMS[1] + CENTNO * CHIMS[2])))

        DXMDT = (CHIMS[1] + 2. * CENTNO * CHIMS[2]) / (86400.0 * CENTCT)

        MANOM0 = self._ZSODSRadto2Pai(MEARTH[0] + self._ZSODSRadto2Pai(CENTNO * (MEARTH[1] + CENTNO * MEARTH[2])))

        DMDT = (MEARTH[1] + 2. * CENTNO * MEARTH[2]) / (86400.0 * CENTCT)

        OMEGA0 = self._ZSODSRadto2Pai(OMEGA[0] + self._ZSODSRadto2Pai(CENTNO * (OMEGA[1] + CENTNO * OMEGA[2])))

        DWDT = (OMEGA[1] + 2. * CENTNO * OMEGA[2]) / (86400.0 * CENTCT)

        EPSLN0 = EPSLON[0] + CENTNO * EPSLON[1]

        DEDT = EPSLON[1] / (86400.0 * CENTCT)

        DELP0 = DELPSI[0] * np.sin(OMEGA0) + DELPSI[1] * np.sin(2.0 * CHIMO)

        DDPDT = DELPSI[0] * np.cos(OMEGA0) * DWDT + DELPSI[1] * np.cos(2. * CHIMO) *2.0 * DXMDT

        DELE0 = DELEPS[0] * np.cos(OMEGA0) + DELEPS[1] * np.cos(2.0 * CHIMO)

        DDEDT = -DELEPS[0] * np.sin(OMEGA0) * DWDT - DELEPS[1] * np.sin(2.0 * CHIMO) * 2.0 * DXMDT

        CHI0 = CHIMO + KAPPA * np.sin(MANOM0)

        DXDT = DXMDT + KAPPA * np.cos(MANOM0) * DMDT

        ANUT0 = DELP0 * np.cos(EPSLN0 + DELE0)

        DANDT = DDPDT * ANUT0 / DELP0 - DELP0 * np.sin(EPSLN0 + DELE0) * (DEDT + DDEDT)

        Dtime = stime - EPTIME

        ALPHA = self._ZSODSRadto2Pai(ALPHA0 + ANUT0 + (DADT + DANDT) * Dtime)

        ALPDOT = DADT + DANDT

        CHISUN = self._ZSODSRadto2Pai(CHI0 + DXDT * Dtime)

        OBLIQE = EPSLN0 + DELE0 + (DEDT + DDEDT) * Dtime

        return ALPHA, ALPDOT, CHISUN, OBLIQE


    def _ZSODSRadto2Pai(self, dRad) :

        dresult = dRad - DPAI*((int)(dRad/DPAI))

        if dRad < 0.0 :
            dresult += DPAI

        return dresult







