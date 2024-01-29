=================================
fy4调用示例
=================================

FY-4数据辐射定标、裁剪
-----------------------------------------

.. code-block:: python

    from fypy.fy4 import fy4scene

    # filename = r'FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20240129041500_20240129041917_4000M_V0001.HDF'
    # filename = r'FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20240129040000_20240129041459_4000M_V0001.HDF'
    filename = r'FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20240129040000_20240129041459_1000M_V0001.HDF'

    mc_pro = fy4scene(filename=filename)
    band1 = mc_pro.Calibration(filename, bandID=2, fillvalue=65535)

    ds = mc_pro.Clip(band1, extent=[70, 20, 135, 55])
    mc_pro.DS2Tiff(r'test.tif', srcDS=ds)
    mc_pro.DS2Netcdf(r'test.nc', 'data', srcDS=ds)
    mc_pro.DS2Hdf(r'test.hdf', 'data', srcDS=ds)

    filename = r'FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20240129040000_20240129041459_4000M_V0001.HDF'

    x = scene.GetGEOData(filename, 'ColumnNumber')
    y = scene.GetGEOData(filename, 'LineNumber')
    NOMSunAzimuth = scene.GetGEOData(filename, 'NOMSunAzimuth')
    NOMSunZenith = scene.GetGEOData(filename, 'NOMSunZenith')
    NOMSatelliteAzimuth = scene.GetGEOData(filename, 'NOMSatelliteAzimuth')
    NOMSatelliteZenith = scene.GetGEOData(filename, 'NOMSatelliteZenith')
    NOMSunGlintAngle = scene.GetGEOData(filename, 'NOMSunGlintAngle')


FY-4数据经纬度、行列号互转
-----------------------------------------

.. code-block:: python

    from fypy.fy4 import fy4scene

    scene = fy4scene(sublon=104.7, resolution=0.04)
    x, y = np.meshgrid(np.arange(scene.colmax), np.arange(scene.rowmax))

    # xy转lon,lat
    lon, lat = scene.xytolatlon(x, y)

    # lat, lon转xy
    x, y = scene.latlon2xy(lat, lon)

FY-4各类角度计算
-----------------------------------------

.. code-block:: python

    from fypy.fy4 import fy4scene

    filename = r'FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20240129040000_20240129041459_4000M_V0001.HDF'
    scene = fy4scene(filename=filename)
    nowdate = scene.StartTime
    print(nowdate)

    lat = np.array([20])
    lon = np.array([104.7])

    suna, sunz = scene.CalSolarAngle(lat, lon, nowdate)
    sata, satz = scene.CalSatAngle(lat, lon)

    # 相对方位角
    rela = scene.CalRelativeAzimuth(sata, suna)

    # 耀斑角
    sungl = scene.CalSunGL(satz, sunz, rela)

