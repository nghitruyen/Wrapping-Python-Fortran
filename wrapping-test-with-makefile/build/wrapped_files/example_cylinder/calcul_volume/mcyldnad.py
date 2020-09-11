"""
Module mcyldnad


Defined at cyldnad.f90 lines 1-10

"""
from __future__ import print_function, absolute_import, division
import _calcul_volume
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

def cyldnad(radius, height):
    """
    vol = cyldnad(radius, height)
    
    
    Defined at cyldnad.f90 lines 4-9
    
    Parameters
    ----------
    radius : Dual_Num
    height : Dual_Num
    
    Returns
    -------
    vol : Dual_Num
    
    """
    vol = _calcul_volume.f90wrap_cyldnad(radius=radius._handle, \
        height=height._handle)
    vol = f90wrap.runtime.lookup_class("calcul_volume.DUAL_NUM").from_handle(vol)
    return vol


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "mcyldnad".')

for func in _dt_array_initialisers:
    func()
