=================================
fy3调用示例
=================================

FY3 MERSI L1数据辐射定标
-----------------------------------------
支持对FY3 MERSI 1KM、250M L1数据辐射定标

.. code-block:: python

    from fypy.fy3 import FY3Orbit, fy3pro
    from fypy.tools import readhdf

    l1file = r'./FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = r'./FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'

    mpro = fy3pro(satID='FY3D', instID='MERSI')
    # 1-4波段
    band14 = mpro.calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')

    # 5 - 19波段
    band519 = mpro.calibration(l1file, '/Data/EV_1KM_RefSB')

    # 20 - 23
    band2023 = mpro.calibration(l1file, '/Data/EV_1KM_Emissive')

    # 24-25
    band2425 = mpro.calibration(l1file, '/Data/EV_250_Aggr.1KM_Emissive')


FY3 轨道数据产品投影
-----------------------------------------
支持对FY3 MERSI L1、L2等轨道产品（ORBT）数据进行投影（WGS84）

.. code-block:: python

    from fypy.fy3 import FY3Orbit, fy3pro
    from fypy.tools import readhdf

    l1file = r'./FY3D_MERSI_GBAL_L1_20220718_0725_1000M_MS.HDF'
    geofile = r'./FY3D_MERSI_GBAL_L1_20220718_0725_GEO1K_MS.HDF'

    mpro = fy3pro(satID='FY3D', instID='MERSI')
    band14 = mpro.calibration(l1file, '/Data/EV_250_Aggr.1KM_RefSB')
    band519 = mpro.calibration(l1file, '/Data/EV_1KM_RefSB')

    lat = readhdf(geofile, '/Geolocation/Latitude')
    lon = readhdf(geofile, '/Geolocation/Longitude')
    mpro = FY3Orbit(satID='FY3D', instID='MERSI')
    dstdata, trans, prj = mpro.translate(data, lat, lon,
                                        resolution=0.01,
                                        vmin=0, vmax=10000,
                                        resampleAlg='near')

    mpro.to_netcdf(r'./test.nc')
    mpro.to_hdf(r'./test.hdf')
    mpro.to_tiff(r'./test.tif')

FY3 10度块产品拼接（支持GLL、HAM投影）
-----------------------------------------
FY3的10度块产品大部分是GLL（等经纬WGS84），由于NDVI产品是HAM投影，新增了对HAM投影

.. code-block:: python

    from fypy.tools import readhdf
    from fypy.fy3.fy3Block10 import FY3Block10

    filelist  = glob.glob(r'./data/NDVI/*.HDF')

    mc_pro = FY3Block10(satID='FY3D', instID='MERSI')
    # 对10度块产品进行拼接
    mc_pro.MergeBlock(filelist, productName='NDVI', sdsname='1000M_10day_NDVI')

    # 保存为TIFF文件
    mc_pro.SaveTiff(r'./test.tif')

    # 将Harmmer投影转换为WGS84
    data, trans, prj = mc_pro.Hammer2Wgs84(mc_pro.dset, resolution=0.01)


