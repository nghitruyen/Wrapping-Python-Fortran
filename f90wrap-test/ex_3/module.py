from __future__ import print_function, absolute_import, division
import _module
import f90wrap.runtime
import logging

class Dual_Num_Auto_Diff(f90wrap.runtime.FortranModule):
    """
    Module dual_num_auto_diff
    
    
    Defined at DNAD.f90 lines 108-1343
    
    """
    @f90wrap.runtime.register_class("module.DUAL_NUM")
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
            result = _module.f90wrap_dual_num_initialise()
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
                _module.f90wrap_dual_num_finalise(this=self._handle)
        
        @property
        def x_ad_(self):
            """
            Element x_ad_ ftype=real(dbl_ad) pytype=float
            
            
            Defined at DNAD.f90 line 119
            
            """
            return _module.f90wrap_dual_num__get__x_ad_(self._handle)
        
        @x_ad_.setter
        def x_ad_(self, x_ad_):
            _module.f90wrap_dual_num__set__x_ad_(self._handle, x_ad_)
        
        @property
        def xp_ad_(self):
            """
            Element xp_ad_ ftype=real(dbl_ad) pytype=float
            
            
            Defined at DNAD.f90 line 120
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _module.f90wrap_dual_num__array__xp_ad_(self._handle)
            if array_handle in self._arrays:
                xp_ad_ = self._arrays[array_handle]
            else:
                xp_ad_ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _module.f90wrap_dual_num__array__xp_ad_)
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
        
    
    @property
    def ndv_ad(self):
        """
        Element ndv_ad ftype=integer(2) pytype=int
        
        
        Defined at DNAD.f90 line 110
        
        """
        return _module.f90wrap_dual_num_auto_diff__get__ndv_ad()
    
    def __str__(self):
        ret = ['<dual_num_auto_diff>{\n']
        ret.append('    ndv_ad : ')
        ret.append(repr(self.ndv_ad))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

dual_num_auto_diff = Dual_Num_Auto_Diff()

class Mcyldnad(f90wrap.runtime.FortranModule):
    """
    Module mcyldnad
    
    
    Defined at cyldnad.f90 lines 1-10
    
    """
    @staticmethod
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
        vol = _module.f90wrap_cyldnad(radius=radius._handle, height=height._handle)
        vol = f90wrap.runtime.lookup_class("module.DUAL_NUM").from_handle(vol)
        return vol
    
    _dt_array_initialisers = []
    

mcyldnad = Mcyldnad()

