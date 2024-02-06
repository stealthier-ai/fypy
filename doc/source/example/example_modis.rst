 =================================
modis调用示例
=================================

modis
-----------------------------------------

.. code-block:: python

    from fypy.modis import modisscene
    from fypy.tools import readhdf4
    
    scene = modisscene()
    filename = r'MOD021KM.A2022001.0330.061.2022001133325.hdf'
    data = readhdf4(filename, 'EV_500_Aggr1km_RefSB')

    geofile = r'MOD03.A2022001.0330.061.2022001090339.hdf'
    lat = readhdf4(geofile, 'Latitude')
    lon = readhdf4(geofile, 'Longitude')

    ds = scene.Reprojection(data, lat, lon,
                            resolution=0.01, resampleAlg='near')
    scene.DS2Tiff(r'test3.tif', ds)

    ds = scene.QuickProject(filename, 'EV_500_Aggr1km_RefSB')
    scene.DS2Tiff(r'test4.tif', ds)

modis L3级数据产品拼接、投影、裁剪
-----------------------------------------

.. code-block:: python

    filename = glob.glob(r'MOD13A2.A2021353.*.hdf')
    sdsname = '1 km 16 days NDVI'
    from fypy.modis import modisscene

    scene = modisscene()
    ds = scene.ConverModisByGDAL(filename, sdsname, resolution=None)
    scene.DS2Tiff(r'test.tif', ds)