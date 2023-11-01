# -*- coding:utf-8 -*-
'''
@Project  : fypy
@File     : __init__.py
@Modify Time      @Author    @Version    
--------------    -------    --------    
2022/7/20 11:29      Lee       1.0         
@Description
------------------------------------
 
'''

from .version import get_versions
__version__ = get_versions()['version']
del get_versions


from fypy.fy3 import fy3L1Pro, FY3Orbit, FY3Block10
from fypy.fy4 import fy4scene
from fypy.draw import colorbar, drawThematic
from fypy.h8 import hsd2hdf, drawH8TrueColor
from fypy.modis import modis2tif

from fypy.RSDP import AtmCorr_FY3D_MERSI, AtmCorr_GF, \
    AtmCorr_Landsat, AtmCorr_Sentinel
