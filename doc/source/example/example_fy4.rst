=================================
fy4调用示例
=================================

fy4pro
-----------------------------------------

.. code-block:: python

    from fypy.fy4 import fy4scene

    geoname = r'./FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20231009041500_20231009041917_2000M_V0001.HDF'

    mc_pro = fy4scene(subpoint=104.7, resolution=0.02)
    satA = mc_pro.getL1Data(geoname, bandID=1)
    # 行列号转经纬度（行列号索引从0开始）
    lat,lon = mc_pro.xytolatlon([2000], [2000])

    # 经纬度转行列号（行列号索引从0开始）
    x, y = mc_pro.latlon2xy(lat, lon)

    mc_pro.clip(satA, dstNodata=65535, outputBounds=[70, 20, 135, 55])
    mc_pro.to_tiff(r'./test.tif')



