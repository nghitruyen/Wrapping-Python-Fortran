"""
Module dual_num_auto_diff


Defined at DNAD.f90 lines 108-1343

"""
from __future__ import print_function, absolute_import, division
import _calcul_volume
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("calcul_volume.DUAL_NUM")
class DUAL_NUM(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=dual_num)
    
    
    Defined at DNAD.f90 lines 115-120
    
    """
    def __init__(self, handle=None):
        """
        self = Dual_Num()
        
        
        Defined at DNAD.f90 lines 115-120
        
        
        Returns
        -------
        this : Dual_Num
        	Object to be constructed
        
        
        Automatically generated constructor for dual_num
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _calcul_volume.f90wrap_dual_num_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Dual_Num
        
        
        Defined at DNAD.f90 lines 115-120
        
        Parameters
        ----------
        this : Dual_Num
        	Object to be destructed
        
        
        Automatically generated destructor for dual_num
        """
        if self._alloc:
            _calcul_volume.f90wrap_dual_num_finalise(this=self._handle)
    
    @property
    def x_ad_(self):
        """
        Element x_ad_ ftype=real(dbl_ad) pytype=float
        
        
        Defined at DNAD.f90 line 119
        
        """
        return _calcul_volume.f90wrap_dual_num__get__x_ad_(self._handle)
    
    @x_ad_.setter
    def x_ad_(self, x_ad_):
        _calcul_volume.f90wrap_dual_num__set__x_ad_(self._handle, x_ad_)
    
    @property
    def xp_ad_(self):
        """
        Element xp_ad_ ftype=real(dbl_ad) pytype=float
        
        
        Defined at DNAD.f90 line 120
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _calcul_volume.f90wrap_dual_num__array__xp_ad_(self._handle)
        if array_handle in self._arrays:
            xp_ad_ = self._arrays[array_handle]
        else:
            xp_ad_ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _calcul_volume.f90wrap_dual_num__array__xp_ad_)
            self._arrays[array_handle] = xp_ad_
        return xp_ad_
    
    @xp_ad_.setter
    def xp_ad_(self, xp_ad_):
        self.xp_ad_[...] = xp_ad_
    
    def __str__(self):
        ret = ['<dual_num>{\n']
        ret.append('    x_ad_ : ')
        ret.append(repr(self.x_ad_))
        ret.append(',\n    xp_ad_ : ')
        ret.append(repr(self.xp_ad_))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def get_ndv_ad():
    """
    Element ndv_ad ftype=integer(2) pytype=int
    
    
    Defined at DNAD.f90 line 110
    
    """
    return _calcul_volume.f90wrap_dual_num_auto_diff__get__ndv_ad()

NDV_AD = get_ndv_ad()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "dual_num_auto_diff".')

for func in _dt_array_initialisers:
    func()
