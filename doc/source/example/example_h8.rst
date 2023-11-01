=================================
h8调用示例
=================================

ahi8_l1_pro
-----------------------------------------

.. code-block:: python

    from fypy.fypy import hsd2hdf
    outdir = r'./H8/20211012'
    hsdpath = r'./H8'
    nowdate = datetime.datetime.strptime('20211012_0320', '%Y%m%d_%H%M')

    hsd2hdf(outdir, nowdate, hsdpath)



