=================================
fy3调用示例
=================================

FY3 MERSI L1数据辐射定标
-----------------------------------------
支持对FY3 MERSI 1KM、250M L1数据辐射定标

.. code-block:: python

    from fypy.fy3 import fy3scene
    from fypy.tools import readhdf

    l1file = r'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = r'FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'

    mc_pro = fy3scene(l1file)
    # 1-4波段
    band14 = mc_pro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')

    # 5 - 19波段
    band519 = mc_pro.Calibration(l1file, '/Data/EV_1KM_RefSB')

    # 20 - 23
    band2023 = mc_pro.Calibration(l1file, '/Data/EV_1KM_Emissive')

    # 24-25
    band2425 = mc_pro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_Emissive')


FY3 轨道数据产品投影
-----------------------------------------
支持对FY3 MERSI L1、L2等轨道产品（ORBT）数据进行投影（WGS84）

.. code-block:: python

    from fypy.fy3 import fy3scene
    from fypy.tools import readhdf

    l1file = r'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = r'FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'
    mc_pro = fy3scene(l1file)
    # 读取通道数据，并进行辐射定标，如果不做辐射定标，使用readhdf方式读取源数据
    data = mc_pro.Calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')

    lat = readhdf(geofile, '/Geolocation/Latitude')
    lon = readhdf(geofile, '/Geolocation/Longitude')

    # 轨道数据投影
    ds = mc_pro.Reprojection(data, lat, lon,
                   resolution=None, vmin=0, vmax=10000, resampleAlg='near')

    mc_pro.DS2Tiff(r'test.tif', ds)

FY3 10度块产品拼接（支持GLL、HAM投影）
-----------------------------------------
FY3的10度块产品大部分是GLL（等经纬WGS84），由于NDVI产品是HAM投影，新增了对HAM投影

.. code-block:: python

    from fypy.fy3 import fy3scene

    filelist  = glob.glob(r'*NDVI*.HDF')

    mc_pro = fy3scene()
    # 对10度块产品进行拼接
    ds = mc_pro.Block10Merge(filelist, ProdID='NDVI', sdsname='1000M_10day_NDVI', fillvalue=-32768)
    mc_pro.DS2Tiff('test1.tif', ds)

    # 将Harmmer投影转换为WGS84
    harmds = mc_pro.Hammer2Wgs84(ds, xRes=0.01, yRes=0.01)
    mc_pro.DS2Tiff('test2.tif', harmds)

FY3 绘制真彩图
-----------------------------------------

.. code-block:: python

    from fypy.fy3 import fy3scene

    l1file = r'FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    scene = fy3scene(filename=l1file)
    scene.load(l1file, ProdID='truecolor')
    scene.show()
    scene.SaveThematic(r'test.png')
