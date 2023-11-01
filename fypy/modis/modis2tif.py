# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : modis2tif.py

@Modify Time : 2022/8/11 15:34

@Author : Lee

@Version : 1.0

@Description :
针对MODIS产品数据，将数据进行投影转换成GeoTiff文件

'''
import glob

# import gdal, osr
import numpy as np
import os
try:
    from osgeo import gdal, gdalconst, osr, ogr
except ImportError:
    try:
        import gdal, gdalconst, osr, ogr
    except ImportError:
        raise ImportError('Python GDAL library not found, please install '
                          'python-gdal')


RESAM_GDAL = ['AVERAGE', 'BILINEAR', 'CUBIC', 'CUBIC_SPLINE', 'LANCZOS',
              'MODE', 'NEAREST_NEIGHBOR']
SINU_WKT = 'PROJCS["Sinusoidal_Sanson_Flamsteed",GEOGCS["GCS_Unknown",' \
           'DATUM["D_unknown",SPHEROID["Unknown",6371007.181,"inf"]],' \
           'PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]' \
           ',PROJECTION["Sinusoidal"],PARAMETER["central_meridian",0],' \
           'PARAMETER["false_easting",0],PARAMETER["false_northing",0]' \
           ',UNIT["Meter",1]]'



class Hdf2TifFor5min():

    def __init__(self, outName, filename, sdsname, resolution=0.01, fillvalue=-999.0):
        datasets = gdal.Open(filename)
        #打开子数据集
        layers = datasets.GetSubDatasets()

        layer = self.GetLayer(layers, sdsname)
        #打开ndvi数据
        raster = gdal.Open(layer)
        #获取元数据
        Metadata = datasets.GetMetadata()
        data = raster.ReadAsArray()

        data = data * 0.001 * 10 + 0.0
        data[data<0] = fillvalue

        #  获取四个角的维度
        Latitudes = Metadata["GRINGPOINTLATITUDE.1"]
        #  采用", "进行分割
        LatitudesList = Latitudes.split(", ")
        #  获取四个角的经度
        Longitude = Metadata["GRINGPOINTLONGITUDE.1"]
        #  采用", "进行分割
        LongitudeList = Longitude.split(", ")

        # 图像四个角的地理坐标
        GeoCoordinates = np.zeros((4, 2), dtype = "float32")
        GeoCoordinates[0] = np.array([float(LongitudeList[0]),float(LatitudesList[0])])
        GeoCoordinates[1] = np.array([float(LongitudeList[1]),float(LatitudesList[1])])
        GeoCoordinates[2] = np.array([float(LongitudeList[2]),float(LatitudesList[2])])
        GeoCoordinates[3] = np.array([float(LongitudeList[3]),float(LatitudesList[3])])

        #  列数
        # Columns = float(Metadata["DATACOLUMNS"])
        Rows, Columns = data.shape
        #  行数
        # Rows = float(Metadata["DATAROWS"])
        #  图像四个角的图像坐标
        PixelCoordinates = np.array([[0, 0],
                                     [Columns - 1, 0],
                                     [Columns - 1, Rows - 1],
                                     [0, Rows - 1]], dtype = "float32")

        #  计算仿射变换矩阵
        from scipy.optimize import leastsq
        def func(i):
            Transform0, Transform1, Transform2, Transform3, Transform4, Transform5 = i[0], i[1], i[2], i[3], i[4], i[5]
            return [Transform0 + PixelCoordinates[0][0] * Transform1 + PixelCoordinates[0][1] * Transform2 - GeoCoordinates[0][0],
                    Transform3 + PixelCoordinates[0][0] * Transform4 + PixelCoordinates[0][1] * Transform5 - GeoCoordinates[0][1],
                    Transform0 + PixelCoordinates[1][0] * Transform1 + PixelCoordinates[1][1] * Transform2 - GeoCoordinates[1][0],
                    Transform3 + PixelCoordinates[1][0] * Transform4 + PixelCoordinates[1][1] * Transform5 - GeoCoordinates[1][1],
                    Transform0 + PixelCoordinates[2][0] * Transform1 + PixelCoordinates[2][1] * Transform2 - GeoCoordinates[2][0],
                    Transform3 + PixelCoordinates[2][0] * Transform4 + PixelCoordinates[2][1] * Transform5 - GeoCoordinates[2][1],
                    Transform0 + PixelCoordinates[3][0] * Transform1 + PixelCoordinates[3][1] * Transform2 - GeoCoordinates[3][0],
                    Transform3 + PixelCoordinates[3][0] * Transform4 + PixelCoordinates[3][1] * Transform5 - GeoCoordinates[3][1]]
        #  最小二乘法求解
        GeoTransform = leastsq(func,np.asarray((1,1,1,1,1,1)))

        #  获取数据时间
        # date = Metadata["RANGEBEGINNINGDATE"]

        #  第一个子数据集合,也就是NDVI数据
        # DatasetNDVI = datasets.GetSubDatasets()[0][0]
        # RasterNDVI = gdal.Open(DatasetNDVI)
        # NDVI = ndviRaster.ReadAsArray()
        # print(tuple(GeoTransform[0]), data.shape)

        self.array2raster(outName, GeoTransform[0], data)
        print(outName,"Saved successfully!")


        #命名输出完整路径文件名

        #进行几何校正
        # geoData = gdal.Warp(outName, ndviRaster,
        #                     dstSRS = 'EPSG:4326', format = 'GTiff',
        #                     resampleAlg = gdal.GRA_Bilinear)
        # del geoData



    def GetLayer(self, layers, sdsname):
        '''
        获取指定的图层的索引名
        :param layers: tuple
        :return: str
        '''

        if sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')
                if sdsname == l_name:
                    return layer[0]

        return None

    #  数组保存为tif
    def array2raster(self, outname, GeoTransform, array, fillvalue=-999.0):
        cols = array.shape[1]  # 矩阵列数
        rows = array.shape[0]  # 矩阵行数
        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(outname, cols, rows, 1, gdal.GDT_Float32)
        # 括号中两个0表示起始像元的行列号从(0,0)开始
        outRaster.SetGeoTransform(tuple(GeoTransform))
        # 获取数据集第一个波段，是从1开始，不是从0开始
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        # 代码4326表示WGS84坐标
        outRasterSRS.ImportFromEPSG(4326)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

        if fillvalue:
            outRaster.GetRasterBand(1).SetNoDataValue(float(fillvalue))


class ConverModisByGDAL():

    def __init__(self, outname, hdfname, sdsname, resolution=None, format="GTiff",
                 epsg=None, wkt=None, resampl='NEAREST_NEIGHBOR', vrt=False):
        """Function for the initialize the object"""
        # Open source dataset

        self.tempfile = None
        self.outname = outname
        self.sdsname = sdsname
        self.resolution = resolution

        # 设置输出文件投影
        if epsg:
            self.dst_srs = osr.SpatialReference()
            self.dst_srs.ImportFromEPSG(int(epsg))
            self.dst_wkt = self.dst_srs.ExportToWkt()
        elif wkt:
            try:
                f = open(wkt)
                self.dst_wkt = f.read()
                f.close()
            except:
                self.dst_wkt = wkt
        else:
            raise Exception('You have to set one of the following option: '
                            '"epsg", "wkt"')
        # error threshold the same value as gdalwarp
        self.maxerror = 0.125
        self.resampling = getResampling(resampl)

        self.driver = gdal.GetDriverByName(format)
        self.vrt = vrt
        if self.driver is None:
            raise Exception('Format driver %s not found, pick a supported '
                            'driver.' % format)

        if isinstance(hdfname, list):
            self._MultiFiles(hdfname)
        elif isinstance(hdfname, str):
            self._GetSourceInfo(hdfname)
            self._createWarped(self.src_driver)
            # self._reprojectOne(self.src_driver)
        else:
            raise Exception('Type for subset parameter not supported')

    def _GetSourceInfo(self, filename, ):
        self.in_name = filename
        self.src_ds = gdal.Open(self.in_name)
        layers = self.src_ds.GetSubDatasets()

        # 获取sdsname所在的图层栅格索引
        self.src_raster = self.GetLayer(layers)
        if self.src_raster is None :
            raise Exception('dataset[%s] is not in the %s' %(self.sdsname, filename))
        self.src_driver = gdal.Open(self.src_raster)
        self.src_proj = self.src_driver.GetProjection()
        self.src_trans = self.src_driver.GetGeoTransform()

        # print(self.src_proj)
        # print(self.src_trans)

        self.src_meta = self.src_driver.GetMetadata()
        #  列数
        # self.tileColumns = int(self.src_meta["DATACOLUMNS"])
        #  行数
        # self.tileRows = int(self.src_meta["DATAROWS"])

        self.tiledata = self.src_driver.ReadAsArray()
        # if 'scale_factor' in self.src_meta and 'add_offset' in self.src_meta :
        #     print('scale_factor: ', np.float(self.src_meta['scale_factor']) , ' add_offset: ', np.float(self.src_meta['add_offset']))
        #     self.tiledata = self.tiledata * np.float32(self.src_meta['scale_factor']) + np.float32(self.src_meta['add_offset'])

        datasize = self.tiledata.shape
        if len(datasize) == 2 :
            #  列数
            self.tileColumns = int(datasize[1])
            #  行数
            self.tileRows = int(datasize[0])
        elif len(datasize) == 3 :
            #  列数
            self.tileColumns = int(datasize[2])
            #  行数
            self.tileRows = int(datasize[1])

        # self.src_meta = src_driver.GetMetadata()
        band = self.src_driver.GetRasterBand(1)

        if '_FillValue' in list(self.src_meta.keys()):
            self.data_fill_value = self.src_meta['_FillValue']
        elif band.GetNoDataValue():
            self.data_fill_value = band.GetNoDataValue()
        else:
            self.data_fill_value = None
        self.datatype = band.DataType


    def _MultiFiles(self, filelist):

        alltrans = []
        alldata = []
        tileindex = []

        countfile = len(filelist)
        if countfile == 0 :
            return None
        elif countfile == 1 :
            self._GetSourceInfo(filelist[0])
            self._createWarped(self.src_driver)
            # self._reprojectOne(self.src_driver)
        else:

            for filename in filelist :
                basename = os.path.basename(filename)
                namelist = basename.split('.')
                if len(namelist[2]) != 6:
                    raise Exception('The row and col id [%s] error ! ' %(namelist[2]))

                tileindex.append([int(namelist[2][1:3]), int(namelist[2][4:]),])

                self._GetSourceInfo(filename)
                alldata.append(self.tiledata)
                alltrans.append(self.src_trans)

                # print(np.mean(self.tiledata), np.max(self.tiledata))
                # print(self.src_trans)

            self.tileindex = np.array(tileindex)
            alldata = np.array(alldata)
            alltrans = np.array(alltrans)

            rowindmax =  np.max(self.tileindex[:, 1])
            rowindmin =  np.min(self.tileindex[:, 1])
            colindmax =  np.max(self.tileindex[:, 0])
            colindmin =  np.min(self.tileindex[:, 0])


            rowtile = int(rowindmax - rowindmin + 1)
            coltile = int(colindmax - colindmin + 1)

            xtotal = int(coltile * self.tileColumns)
            ytotal = int(rowtile * self.tileRows)

            dtype = self._GetNumpyType(self.datatype)
            self.srcdata = np.full(shape=(ytotal, xtotal), fill_value=self.data_fill_value, dtype=dtype)

            for i in range(countfile) :
                indx = self.tileindex[i][0]
                indy = self.tileindex[i][1]
                offx = (indx - colindmin) * self.tileColumns
                offy = (indy - rowindmin) * self.tileRows

                self.srcdata[offy:offy+self.tileRows, offx:offx+self.tileColumns] = alldata[i]

            self.src_trans = tuple([np.min(alltrans[:, 0]), np.mean(alltrans[:, 1]), np.mean(alltrans[:, 2]),
                                   np.max(alltrans[:, 3]), np.mean(alltrans[:, 4]), np.mean(alltrans[:, 5])])

            name = os.path.basename(self.outname)
            pathdir = os.path.dirname(self.outname)
            self.tempfile = os.path.join(pathdir, str(hash(name)) + '.tif')
            dr = gdal.GetDriverByName("GTiff")
            self.src_driver = dr.Create(self.tempfile, xtotal, ytotal, 1, self._GetGdalType(dtype))

            self.src_driver.SetProjection(self.src_proj)
            self.src_driver.SetGeoTransform(self.src_trans)
            self.src_driver.GetRasterBand(1).WriteArray(self.srcdata)
            # if self.data_fill_value:
            #     self.src_driver.GetRasterBand(1).SetNoDataValue(float(self.data_fill_value))
            #     self.src_driver.GetRasterBand(1).Fill(float(self.data_fill_value))

            self._createWarped(self.src_driver)

    def __del__(self):
        if self.tempfile :
            tempfile = self.tempfile

            del self.src_driver
            if os.path.isfile(tempfile) :
                os.remove(tempfile)
                tempfile = None


    def _createWarped(self, src_driver):
        '''
        Create a warped VRT file to fetch default values for target raster
        dimensions and geotransform

        :param raster: the name of raster, for HDF have to be one subset
        :return:
        '''

        # 投影转换
        tmp_ds = gdal.AutoCreateWarpedVRT(src_driver, self.src_proj,
                                          self.dst_wkt, self.resampling,
                                          self.maxerror)

        if not self.resolution:
            self.dst_xsize = tmp_ds.RasterXSize
            self.dst_ysize = tmp_ds.RasterYSize
            self.dst_trans = tmp_ds.GetGeoTransform()
        else:
            bbox = self._boundingBox(tmp_ds)
            self.dst_xsize = self._calculateRes(bbox[0][0], bbox[1][0],
                                                self.resolution)
            self.dst_ysize = self._calculateRes(bbox[0][1], bbox[1][1],
                                                self.resolution)
            if self.dst_xsize == 0:
                raise Exception('Invalid number of pixel 0 for X size. The '
                                'problem could be in an invalid value of '
                                'resolution')
            elif self.dst_ysize == 0:
                raise Exception('Invalid number of pixel 0 for Y size. The '
                                'problem could be in an invalid value of '
                                'resolution')
            self.dst_trans = [bbox[0][0], self.resolution, 0.0,
                              bbox[1][1], 0.0, -self.resolution]
        del tmp_ds


    def _boundingBox(self, src):
        """Obtain the bounding box of raster in the new coordinate system

           :param src: a GDAL dataset object

           :return: a bounding box value in lists
        """
        src_gtrn = src.GetGeoTransform(can_return_null=True)

        src_bbox_cells = ((0., 0.),
                          (0, src.RasterYSize),
                          (src.RasterXSize, 0),
                          (src.RasterXSize, src.RasterYSize))

        geo_pts_x = []
        geo_pts_y = []
        for x, y in src_bbox_cells:
            x2 = src_gtrn[0] + src_gtrn[1] * x + src_gtrn[2] * y
            y2 = src_gtrn[3] + src_gtrn[4] * x + src_gtrn[5] * y
            geo_pts_x.append(x2)
            geo_pts_y.append(y2)
        return ((min(geo_pts_x), min(geo_pts_y)), (max(geo_pts_x),
                                                   max(geo_pts_y)))


    def _calculateRes(self, minn, maxx, res):
        """Calculate the number of pixel from extent and resolution

           :param float minn: minimum value of extent
           :param float maxx: maximum value of extent
           :param int res: resolution of output raster

           :return: integer number with the number of pixels
        """
        return int(round((maxx - minn) / res))


    def _progressCallback(self, pct, message, user_data):
        """For the progress status"""
        return 1  # 1 to continue, 0 to stop


    def _GetGdalType(self, datatype):
        if datatype == np.byte :
            return gdal.GDT_Byte
        elif datatype == np.uint16 :
            return gdal.GDT_UInt16
        elif datatype == np.int16 :
            return gdal.GDT_Int16
        elif datatype == np.uint32 :
            return gdal.GDT_UInt32
        elif datatype == np.int32 :
            return gdal.GDT_Int32
        elif datatype == np.float32 :
            return gdal.GDT_Float32
        elif datatype == np.float64 :
            return gdal.GDT_Float64
        else:
            return gdal.GDT_Unknown

    def _GetNumpyType(self, datatype):
        if datatype == gdal.GDT_Byte :
            return np.byte
        elif datatype == gdal.GDT_UInt16 :
            return np.uint16
        elif datatype == gdal.GDT_Int16 :
            return np.int16
        elif datatype == gdal.GDT_UInt32 :
            return np.uint32
        elif datatype == gdal.GDT_Int32 :
            return np.int32
        elif datatype == gdal.GDT_Float32 :
            return np.float32
        elif datatype == gdal.GDT_Float64 :
            return np.float64
        else:
            return None

    def GetLayer(self, layers):
        '''
        获取指定的图层的索引名
        :param layers: tuple
        :return: str
        '''

        if self.sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')
                # print(self.sdsname, l_name)
                if self.sdsname == l_name:
                    return layer[0]

        return None

    def savetif(self, outname):

        out_name = outname

        if self.vrt:
            out_name = "{pref}.tif".format(pref=self.outname)

        # 创建Tiff文件
        try:
            dst_ds = self.driver.Create(out_name,
                                        self.dst_xsize, self.dst_ysize,
                                        1, self.datatype)
        except:
            raise Exception('Not possible to create dataset %s' % out_name)

        dst_ds.SetProjection(self.dst_wkt)
        dst_ds.SetGeoTransform(self.dst_trans)

        if self.data_fill_value:
            dst_ds.GetRasterBand(1).SetNoDataValue(float(self.data_fill_value))
            dst_ds.GetRasterBand(1).Fill(float(self.data_fill_value))
        cbk = self._progressCallback


        # value for last parameter of above self._progressCallback
        cbk_user_data = None
        try:
            gdal.ReprojectImage(self.src_driver, dst_ds, self.src_proj,
                                self.dst_wkt, self.resampling, 0,
                                self.maxerror, cbk, cbk_user_data)
            # if not quiet:
            #     print("Layer {name} reprojected".format(name=l))
        except:
            raise Exception('Not possible to reproject dataset '
                            '{name}'.format(name=self.sdsname))
        dst_ds.SetMetadata(self.src_meta)

        del dst_ds

        print('save %s success...' %(out_name))


def modis2tif(outname, filename, sdsname, resolution=None,
            format='GTiff', epsg=4326, wkt=None, resampl='NEAREST_NEIGHBOR') :
    '''
    将MODIS产品hdf文件进行投影转换，输出geotif
    可处理MODIS 5分钟段产品和sin grid产品文件
    :param outname: str, 输出geotiff文件名
    :param filename: str or list, 输入待拼接、转换的文件
    :param sdsname: 数据集名
    :param resolution: float，degree， 输出结果的分辨率
    :param format:
    :param epsg:
    :param wkt:
    :param resampl: 重采样方式
    :return:
    '''
    try:
        mds = ConverModisByGDAL(outname, filename, sdsname, resolution=resolution,
                                format=format, epsg=epsg, wkt=wkt, resampl=resampl)
        mds.savetif(outname)
    except BaseException as e :
        print(e)
        return False

    return True


def getResampling(res):
    """Return the GDAL resampling method

       :param str res: the string of resampling method
    """
    if res == 'AVERAGE':
        return gdal.GRA_Average
    elif res == 'BILINEAR' or res == 'BICUBIC':
        return gdal.GRA_Bilinear
    elif res == 'LANCZOS':
        return gdal.GRA_Lanczos
    elif res == 'MODE':
        return gdal.GRA_Mode
    elif res == 'NEAREST_NEIGHBOR':
        return gdal.GRA_NearestNeighbour
    elif res == 'CUBIC_CONVOLUTION' or res == 'CUBIC':
        return gdal.GRA_Cubic
    elif res == 'CUBIC_SPLINE':
        return gdal.GRA_CubicSpline








'''
由于L1A、L1B、Geolocation数据均使用HDF－EOS格式，
因此在介绍L1A、L1B、Geolocation数据格式之前
有必要对HDF－EOS进行简要的说明。
HDF-EOS是NASA为遥感应用而对
NCSA（National Center for Supercomputing Applications
 美国国家超级计算中心）的HDF（Hierarchical Data Format 分级数据格式）进行的扩充。

5.MODIS数据产品分级

5.1 MODIS数据产品分级系统
    MODIS标准数据产品分级系统由5级数据构成，它们分别是：0级、1级、2级、3级和4级

5.2 0级数据
    卫星地面站直接接收到的、未经处理的、包括全部数据信息在内的原始数据为0级数据。

5.3 1级数据
    对没有经过处理的、完全分辨率的仪器数据进行重建，数据时间配准，使用辅助数据注解，
    计算和增补到0级数据之后为1级数据。

5.4 2级数据
    在1级数据基础上开发出的、具有相同空间分辨率和覆盖相同地理区域的数据为2级数据。

5.5 3级数据
    3级数据时以统一的时间-空间栅格表达的变量，通常具有一定的完整性和一致性。
    在3级水平上，将可以集中进行科学研究，如：定点时间序列，来自单一技术的
    观测方程和通用模型等。

5.6 4级数据
    通过分析模型和综合分析3级以下数据得出的结果数据为4级数据。

6 标准数据产品类型
    该部分详细链接https://wenku.baidu.com/view/ffb3c468a98271fe910ef9bf.html

6.1
MODIS 标准数据产品根据内容的不同分为 0 级、1 级数据产品，
在 1B 级数据产品之后，划分 2 －4 级数据产品，
包括：陆地标准数据产品、大气标准数据产品和海洋标准数据产品等三种主要标准数据产品类型，
总计分解为 44 种标准数据产品类型。它们分别是：

MOD01：即 MODIS1A 数据产品。
MOD02：即 MODIS1B 数据产品。
MOD03：即 MODIS 数据地理定位文件。
MOD04：大气 2、3 级标准数据产品，内容为气溶胶产品，
    Lambert 投影空间分辨率 1 公里，地理坐标 30 秒空间分辨率，
    每日数据为 2 级数据产品，每旬、每月数据合成为 3 级数据产品。
MOD05：可降水量。2 级大气产品。
MOD06：大气 2、3 级标准数据产品，内容为云产品，
    Lambert 投影空间分辨率 1 公里，地理坐标 30 秒空间分辨率，
    每日数据为 2 级数据产品，每旬、每月数据合成为 3 级数据产品。
MOD07：大气 2、3 级标准数据产品，内容为大气剖面数据，
    Lambert 投影空间分辨率 1 公里，地理坐 标 30 秒空间分辨率，
    每日数据为 2 级数据产品，每旬、每月数据合成为 3 级数据产品。
MOD08：大气 3 级标准数据产品，内容为栅格大气产品，
    1 公里空间分辨率。每日、每旬、每月合成 数据。
MOD09：陆地 2 级标准数据产品，内容为表面反射；
    空间分辨率 250m；白天每日数据。
MOD10：陆地 2、3 级标准数据产品，内容为雪覆盖，
    每日数据为 2 级数据，空间分辨率 500 米，
    旬、 月数据合成为 3 级数据，空间分辨率 500 米。
MOD11：陆地 2、3 级标准数据产品，内容为地表温度和辐射率，
    Lambert 投影，空间分辨率 1 公里， 
    地理坐标为 30 秒，每日数据为 2 级数据，
    每旬、每月数据合成为 3 级数据。
MOD12：陆地 3 级标准数据产品，
    内容为土地覆盖/土地覆盖变化，1km，1/4?，季节的，
    生物地球化 学循环，土地覆盖变化，3 级数据产品。
MOD13：陆地 2 级标准数据产品，内容为栅格的归一化植被指数和增强型植被指数（NDVI/EVI），空间分辨率 250m。
MOD14：陆地 2 级标准数据产品，内容为热异常-火灾和生物量燃烧，
    空间分辨率 1km，确定火灾发生 的位置、火灾等级以及暗火与燃烧比。
MOD15：陆地 3 级标准数据产品，内容为叶面积指数和光合有效辐射，
    空间分辨率 1km，每天的及旬、 月合成产品。
MOD16：陆地 4 级标准数据产品，内容为蒸腾作用，空间分辨率 1km，旬、月合成产品。
MOD17：陆地 4 级标准数据产品，内容为植被产品，NPP，空间分辨率为 250 米，1 公里，旬、月度 频率。
MOD18：海洋 2、3 级标准数据产品，内容为标准的水面辐射，全球洋面，空间分辨率 1km，日、旬、 月，海洋叶绿素。
MOD19：海洋 2、3 级标准数据产品，内容为色素浓度，全球洋面，空间分辨率 1km，日、旬、月度 数据。
MOD20：海洋 2、3 级标准数据产品，内容为叶绿素荧光性，全球洋面，空间分辨率 
        1km，叶绿素水 平大于 2.0mg/m3，日、旬、月度数据。 第 2 页 共 20 页 MODIS 数据格式
MOD21：海洋 2 级标准数据产品，内容为叶绿素-色素浓度，空间分辨率 1km，日、旬、月度数据。
MOD22：海洋 2、3 级标准数据产品，内容为光合可利用辐射（PAR） ，全球洋面，1km，日、旬、月 度数据。
MOD23：海洋 3 级标准数据产品，内容为悬浮物浓度。
MOD24：海洋 3 级标准数据产品，内容为有机质浓度。
MOD25：海洋 2、3 级标准数据产品，内容为球石浓度，全球洋面，空间分辨率 1km、20km，日、旬、 月度数据。
MOD26：海洋 3 级标准数据产品，内容为海洋水衰减系数。
MOD27：海洋 2、3 级标准数据产品，内容为海洋初级生产力，全球洋面，空间分辨率 1km，日、旬、 月度数据。
MOD28：海洋 2、3 级标准数据产品，内容为海面温度，全球洋面，空间分辨率 1km，每天的，每周 的/昼夜的，能量和水平衡，气候变化模型。
MOD29：海洋 2 级标准数据产品，内容为海冰覆盖，海洋，1 公里分辨率，日、旬数据。
MOD30： （未定）
MOD31：海洋 2、3 级标准数据产品，内容为藻红蛋白浓度， 1 公里分辨率，日、旬、月度数据。
MOD32：海洋 2 级标准数据产品，内容为处理框架和匹配的数据库，1 公里分辨率，日、旬、月度数 据，用于海洋叶绿素、海洋生产力计算。
MOD33：陆地 3 级标准数据产品，内容为雪覆盖，空间分辨率 500 米，日、旬、月度数据。
MOD34： （未定）
MOD35：大气 2 级标准数据产品，内容为云掩膜，空间分辨率为 250m 和 1 公里，日数据。
MOD36：海洋 3 级标准数据产品，内容为总吸收系数，空间分辨率为 1 公里，日、旬、月度数据
MOD37：海洋 2、3 级标准数据产品，内容为海洋气溶胶特性，空间分辨率 1km，日、旬、月度数据。
MOD38： （未定）
MOD39：海洋 2、3 级标准数据产品，内容为纯水势，空间分辨率 1km，日、旬、月度数据。
MOD40：陆地 3 级标准数据产品，内容为栅格的热异常，空间分辨率 1 公里，日、旬、月度数据。
MOD41： （未定）
MOD42：海洋 3 级标准数据产品，内容为海冰覆盖，空间分辨率 1 公里，日、旬、月度数据。
MOD43：陆地 3 级标准数据产品，内容为表面反射，BRDF/Albedo 参数，空间分辨率 1 公里，日、旬、 月度数据。
MOD44：陆地 3 级标准数据产品，内容为植被覆盖转换，250m，季度、年度，判定植被覆盖转换的发 生和类型。

6.2 特殊数据产品类型
MOD45：在 MOD02（1B）数据基础上，经过 BOWTIE 处理后的数据产品。
MOD46：在 MOD02（1B）数据基础上，经过 BOWTIE 处理后，并经过除云后的数据产品。 （其他特 殊数据产品待列） 。

7.MODIS产品MOD13Q1简介
MODIS产品有44种，可以分为大气、陆地、冰雪、海洋四个专题数据产品，
其中MOD13Q1属于陆地专题的产品，全称为MODIS/Terra Vegetation 
Indices 16-Day L3 Global 250m SIN Grid.，简称：MOD13Q1。

全球的MOD13Q1数据是一个采用Sinusoidal投影方式的3级网格数据产品，
具有250米的空间分辨率，每隔16天提供一次。
当缺少250米分辨率的蓝波段时，evi算法使用500米分辨率的蓝波段矫正残余的大气影响。
'''



