# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4core.py

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
            raise Exception("resolution error! please input [0.0025, 0.005, 0.01, 0.02, 0.04]")


        self.FLAT = 0.00335281317789691
        self.AE = 6378.137
        self.E2 = 0.0066943800699785


    def xy2latlon(self, x, y, fillvalue=65535):
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

        lon[flag] = fillvalue
        lat[flag] = fillvalue

        lon[fillflag] = fillvalue
        lat[fillflag] = fillvalue

        lat = np.array(lat, dtype=np.float32)
        lon = np.array(lon, dtype=np.float32)

        return lon, lat

    def latlon2xy(self, lat, lon, fillvalue=65535):
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

        lat1 = np.array(lat.copy())
        lon1 = np.array(lon.copy())

        # lon1[lon1 > 180] -= 360.0
        fillflag = (lat1<-90) | (lat1 > 90) | \
                   (lon1<-180) | (lon1 > 180)

        lon1 = lon1 * deg2rad
        lat1 = lat1 * deg2rad

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

        line[(line<0) | (line >= self.rowmax)] = fillvalue
        col[(col<0) | (col>=self.colmax)]  = fillvalue
        col[fillflag]  = fillvalue
        line[fillflag]  = fillvalue

        return col, line

    def CalSolarAngle(self, lat, lon, nowdate, fillvalue=65535):
        ''' 计算太阳天顶角SOZ和方位角SOA '''

        dHours, dMinutes, dSeconds = nowdate.hour, nowdate.minute, nowdate.second
        iYear, iMonth, iDay = nowdate.year, nowdate.month, nowdate.day

        dEarthMeanRadius = 6371.01
        dAstronomicalUnit = 149597890

        dDecimalHours = dHours + (dMinutes + dSeconds / 60.) / 60.

        liAux1 = int((iMonth - 14.) / 12.)
        liAux2 = int((1461. * (iYear + 4800. + liAux1)) / 4.) + int((367. * (iMonth - 2. - 12. * liAux1)) / 12.) - int(
            (3. * int((iYear + 4900. + liAux1) / 100.)) / 4.) + iDay - 32075.
        dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.

        dElapsedJulianDays = dJulianDate - 2451545.0

        dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
        dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
        dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
        dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) + 0.00034894 * np.sin(
            2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
        dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * np.cos(dOmega)

        dSin_EclipticLongitude = np.sin(dEclipticLongitude)
        dY = np.cos(dEclipticObliquity) * dSin_EclipticLongitude
        dX = np.cos(dEclipticLongitude)
        dRightAscension = np.arctan2(dY, dX)
        if dRightAscension < 0.0:
            dRightAscension = dRightAscension + 2.0 * np.pi
        dDeclination = np.arcsin(np.sin(dEclipticObliquity) * dSin_EclipticLongitude)

        dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours
        dLocalMeanSiderealTime = np.deg2rad(dGreenwichMeanSiderealTime * 15. + lon)
        dHourAngle = dLocalMeanSiderealTime - dRightAscension
        dLatitudeInRadians = np.deg2rad(lat)
        dCos_Latitude = np.cos(dLatitudeInRadians)
        dSin_Latitude = np.sin(dLatitudeInRadians)
        dCos_HourAngle = np.cos(dHourAngle)
        Zenith = (np.arccos(dCos_Latitude * dCos_HourAngle * np.cos(dDeclination)
                                  + np.sin(dDeclination) * dSin_Latitude))
        dY = -np.sin(dHourAngle)
        dX = np.tan(dDeclination) * dCos_Latitude - dSin_Latitude * dCos_HourAngle
        Azimuth = np.arctan2(dY, dX)
        Azimuth[Azimuth < 0.0] = Azimuth[Azimuth < 0.0] + 2.0 * np.pi
        Azimuth = np.rad2deg(Azimuth)
        # Parallax Correction
        dParallax = (dEarthMeanRadius / dAstronomicalUnit) * np.sin(Zenith)
        Zenith = np.rad2deg(Zenith + dParallax)

        flag = (lat <-90) | (lat > 90) |\
               (lon < -180) | (lon > 180)
        Azimuth[flag] = fillvalue
        Zenith[flag] = fillvalue

        return Azimuth, Zenith

    def CalSatAngle(self, lat, lon, fillvalue=65535):
        ''' 计算卫星的天顶角和方位角 '''

        S = self.subpoint / 180. * np.pi
        N = lon / 180. * np.pi
        L = lat / 180. * np.pi
        G = S - N

        LL = np.sqrt(1 - (np.cos(G) * np.cos(L)) ** 2)
        E = np.arctan((np.cos(G) * np.cos(L) - 0.15127) / LL)
        # Azimuth
        A = np.pi - np.arctan2(np.tan(G), np.sin(L))
        # 弧度->度
        Azimuth = A / np.pi * 180
        Elevat = E / np.pi * 180
        Zenith = (90 - Elevat)

        flag = (lat <-90) | (lat > 90) | \
               (lon < -180) | (lon > 180)
        Azimuth[flag] = fillvalue
        Zenith[flag] = fillvalue

        return Azimuth, Zenith

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

    def CalRelativeAzimuth(self, sata, suna, fillvalue=65535):
        fillflag = (sata > 360) | (sata < 0) | (suna > 360) | (suna < 0)

        RelativeAzi = np.fabs(sata - suna) + 180.0 # ???? + 180

        RelativeAzi[fillflag] = fillvalue

        return RelativeAzi


def GetNameInfo(filename):
    ''' 根据输入数据的文件名获取信息 '''
    # FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20211211060000_20211211061459_4000M_V0001.HDF
    # FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20231231040000_20231231041459_4000M_V0001.HDF
    # FY4A-_AGRI--_N_DISK_1047E_L2-_CLM-_MULT_NOM_20240106030000_20240106031459_4000M_V0001.NC
    # FY4A-_AGRI--_N_REGC_1047E_L2-_CLM-_MULT_NOM_20231230011500_20231230011917_4000M_V0001.NC

    nameinfo = {}
    basename = os.path.basename(filename)
    if len(basename) != 88 and len(basename) != 89 :
        print('输入文件名为非标准文件名，将不作信息提取【%s】' %(basename))
        return nameinfo

    basename = basename.replace('-', '')
    namelist = basename.split('_')
    nameinfo['SatID']  = namelist[0]
    nameinfo['InstID'] = namelist[1]
    nameinfo['ObsType'] = namelist[2]
    nameinfo['RegionID'] = namelist[3]
    nameinfo['SubLon'] = namelist[4]
    nameinfo['LevelID'] = namelist[5]
    nameinfo['ProdID'] = namelist[6]
    nameinfo['Proj'] = namelist[8]
    nameinfo['StartTime'] = datetime.datetime.strptime(namelist[9], '%Y%m%d%H%M%S')
    nameinfo['EndTime'] = datetime.datetime.strptime(namelist[10], '%Y%m%d%H%M%S')
    if 'KM' in namelist[11] :
        nameinfo['Resolution'] = float(namelist[11].replace('KM',''))/100.0
    elif 'M' in namelist[11] :
        nameinfo['Resolution'] = float(namelist[11].replace('M',''))/100.0/1000.0
    nameinfo['Version'] = namelist[12]

    return nameinfo




