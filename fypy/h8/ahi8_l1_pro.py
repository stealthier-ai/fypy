# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : ahi8_l1_pro.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee

@Version     : 1.0

@Description :

'''

import os
import datetime
import shutil
import time

import numpy as np
import glob
from osgeo import gdal, osr
import warnings
warnings.filterwarnings('ignore')

from .ahi8_read_hsd import ahi8_read_hsd
from fypy.tools import writehdf


class ahi8_l1_pro(ahi8_read_hsd):

    def __init__(self, filelist, bandid=None, obstype=None, tmppath=None, fillvalue=65535):

        self.BandID = bandid
        self.ObsType = obstype
        self.tmpfile = []
        outdata = None

        for hsdname in filelist :
            if not os.path.isfile(hsdname):
                continue

            print('正在解压bz2文件【%s】' %(hsdname))
            self._unzipped = self.unzip_file(hsdname, tmppath)
            if self._unzipped:
                # But if it is, set the filename to point to unzipped temp file
                self.is_zipped = True
                filename = self._unzipped

                self.tmpfile.append(filename)
            else:
                filename = hsdname

            if filename.endswith('.bz2') :
                print('解压bz2文件失败【%s】' %(filename))
                continue

            # print(filename)

            # 根据块号对数据进行拼接
            SegmentNum, data = self.readhsd(filename)
            if outdata is None :
                line, pixel = data.shape
                outdata = np.full(shape=(line*self.SegmentTotal, pixel), fill_value=fillvalue, dtype=np.uint16)

            data[np.isnan(data)] = fillvalue/100.0
            outdata[(SegmentNum-1)*line:(SegmentNum)*line, :] = np.array(data*100, dtype=np.uint16)

        if outdata is None :
            self.data = None
        else:
            self.data = outdata

    def toGeoTiff(self, data, outname=None, bbox=None, subpoint=140.7,
               srcRes=0.04, dstRes=None,
               fillvalue=None, dstSRS='EPSG:4326'):
        '''
        将DISK 投影转换成等经纬投影

        Parameters
        ----------
        data : numpy.array
            输入数据
        outname: str, optional
            输出文件名
        bbox : tuple, optional
            output bounds as (minX, minY, maxX, maxY) in target SRS
        subpoint: float
            星下点
        resolution: float
            数据图像的分辨率，单位为degree
        fillvalue : float
            数据填充值
        dstSRS

        Returns
        -------

        '''
        if dstRes is None :
            dstRes = srcRes

        im_data = np.array(data, dtype=np.float32)
        dtype = self._gettype(im_data.dtype)
        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        elif len(im_data.shape) == 2:
            im_bands, (im_height, im_width) = 1,im_data.shape
        else:
            im_bands, (im_height, im_width) = 1,im_data.shape

        Driver = gdal.GetDriverByName('MEM')
        memDs = Driver.Create('', im_height, im_width, im_bands, dtype)

        # 写入数据
        if im_bands == 1:
            memDs.GetRasterBand(1).WriteArray(im_data)
        else:
            for i in range(im_bands):
                memDs.GetRasterBand(i+1).WriteArray(im_data[i])

        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 +lon_0={subpoint} +no_defs'.format(
            subpoint=subpoint ))
        memDs.SetProjection(srs.ExportToWkt())
        memDs.SetGeoTransform([-5496000, srcRes*100*1000, 0,
                               5496000, 0, -srcRes*100*1000])
        if outname is None :
            warpDs = gdal.Warp('', memDs, format='MEM', dstSRS=dstSRS,
                               outputBounds=bbox, xRes=dstRes, yRes=dstRes,
                               dstNodata=fillvalue)
        else:
            warpDs = gdal.Warp(outname, memDs, format='GTiff', dstSRS=dstSRS,
                               outputBounds=bbox, xRes=dstRes, yRes=dstRes,
                               srcNodata=fillvalue, dstNodata=fillvalue,
                               creationOptions=["COMPRESS=LZW"])
        if warpDs is None :
            return None

        data_GLL = warpDs.ReadAsArray()
        del warpDs

        return data_GLL

    def readhsd(self, filename):
        ''' 解析hdf文件 '''

        name = os.path.basename(filename)
        namelist = name.split('_')

        # self.NowTime = datetime.datetime.strptime('%s %s' % (namelist[2], namelist[3]), '%Y%m%d %H%M')
        # self.BandID = int(namelist[4][1:])
        # self.ObsType = namelist[5]
        self.Resolution = float(namelist[6][1:])/1000.0
        SegmentNum = int(namelist[7][1:3])
        self.SegmentTotal = int(namelist[7][3:5])

        if self.GetData(filename, SegmentNum, self.SegmentTotal) :
            if self.data is not None :
                data = self.data.values
                return SegmentNum, data
            else:
                return SegmentNum, None
        else:
            return SegmentNum, None

    def unzip_file(self, filename, tmppath=None):
        from subprocess import Popen, PIPE
        from io import BytesIO
        from contextlib import closing
        import shutil
        import bz2

        try:
            from shutil import which
        except ImportError:
            # python 2 - won't be used, but needed for mocking in tests
            which = None

        """Unzip the file if file is bzipped = ending with 'bz2'."""
        try:
            if filename.endswith('bz2'):
                # fdn, tmpfilepath = tempfile.mkstemp()
                if tmppath is None :
                    tmpfilepath = filename.replace('.bz2','')
                else:
                    tmpfilepath = os.path.join(tmppath, os.path.basename(filename).replace('.bz2',''))

                if os.path.isfile(tmpfilepath) :
                    return tmpfilepath
                # print("解压bz2文件【%s】" %(tmpfilepath))
                # try pbzip2
                pbzip = which('pbzip2')
                # Run external pbzip2
                if pbzip is not None:
                    n_thr = os.environ.get('OMP_NUM_THREADS')
                    if n_thr:
                        runner = [pbzip,
                                  '-dc',
                                  '-p'+str(n_thr),
                                  filename]
                    else:
                        runner = [pbzip,
                                  '-dc',
                                  filename]
                    p = Popen(runner, stdout=PIPE, stderr=PIPE)
                    stdout = BytesIO(p.communicate()[0])
                    status = p.returncode
                    if status != 0:
                        raise IOError("pbzip2 error '%s', failed, status=%d"
                                      % (filename, status))
                    with closing(open(tmpfilepath, 'wb')) as ofpt:
                        try:
                            stdout.seek(0)
                            shutil.copyfileobj(stdout, ofpt)
                        except IOError:
                            import traceback
                            traceback.print_exc()
                            print("Failed to read bzipped file %s",
                                  str(filename))
                            os.remove(tmpfilepath)

                    return tmpfilepath

                # Otherwise, fall back to the original method
                bz2file = bz2.BZ2File(filename)
                with closing(open(tmpfilepath, 'wb')) as ofpt:
                    try:
                        ofpt.write(bz2file.read())
                    except IOError:
                        import traceback
                        traceback.print_exc()
                        print("Failed to read bzipped file %s", str(filename))
                        # os.remove(tmpfilepath)
                        return None
                return tmpfilepath
        except BaseException :
            return None

    def __del__(self):
        for filename in self.tmpfile :
            if os.path.isfile(filename) :
                try:
                    os.remove(filename)
                except BaseException as e :

                    time.sleep(2)
                    try:
                        fp = open(filename, 'r')
                        fp.close()
                        os.remove(filename)
                    except BaseException as e :
                        print(e)


def hsd2hdf(outdir, hsdpath, nowdate, SatID='H08', fillvalue=65535):

    outname = os.path.join(outdir, 'HS_H08_%s_FLDK.h5' %(nowdate.strftime('%Y%m%d_%H%M')))
    overwrite = 1

    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    tempdir = os.path.join(outdir, 'temp')
    if not os.path.isdir(tempdir) :
        os.makedirs(tempdir)

    tempfile = []
    for chid in np.arange(1, 17) :
        searchfile = os.path.join(hsdpath,
                                  r'*HS_{satid}_{date}_B{chid}_FLDK_R*_S*.DAT.bz2'.format(
                                      satid=SatID,
                                      date=nowdate.strftime('%Y%m%d_%H%M'),
                                      chid='%02d' %(chid)))
        filelist = glob.glob(searchfile)
        if len(filelist) == 0:
            print('未匹配到文件【%s】' %(searchfile))
            continue

        mpro = ahi8_l1_pro(filelist, bandid=chid, obstype='FLDK', tmppath=tempdir, fillvalue=fillvalue)
        data = mpro.data
        tempfile.extend(mpro.tmpfile)
        del mpro
        mpro=None

        writehdf(outname, 'band%02d' %(chid), data, overwrite=overwrite)
        overwrite = 0

    deletetempfile(tempfile)


def deletetempfile(filelist):
    for filename in filelist :
        if os.path.isfile(filename) :
            try:
                os.remove(filename)
            except BaseException as e :
                time.sleep(2)
                try:
                    fp = open(filename, 'r')
                    fp.close()
                    os.remove(filename)
                except BaseException as e :
                    print(e)