from __future__ import print_function, absolute_import, division
import _module
import f90wrap.runtime
import logging

class Module_Test(f90wrap.runtime.FortranModule):
    """
    Module module_test
    
    
    Defined at ./test.f90 lines 1-11
    
    """
    @f90wrap.runtime.register_class("module.real_array")
    class real_array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=real_array)
        
        
        Defined at ./test.f90 lines 2-3
        
        """
        def __init__(self, handle=None):
            """
            self = Real_Array()
            
            
            Defined at ./test.f90 lines 2-3
            
            
            Returns
            -------
            this : float
            	Object to be constructed
            
            
            Automatically generated constructor for real_array
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _module.f90wrap_real_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Real_Array
            
            
            Defined at ./test.f90 lines 2-3
            
            Parameters
            ----------
            this : float
            	Object to be destructed
            
            
            Automatically generated destructor for real_array
            """
            if self._alloc:
                _module.f90wrap_real_array_finalise(this=self._handle)
        
        @property
        def item(self):
            """
            Element item ftype=real pytype=float
            
            
            Defined at ./test.f90 line 3
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _module.f90wrap_real_array__array__item(self._handle)
            if array_handle in self._arrays:
                item = self._arrays[array_handle]
            else:
                item = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _module.f90wrap_real_array__array__item)
                self._arrays[array_handle] = item
            return item
        
        @item.setter
        def item(self, item):
            self.item[...] = item
        
        def __str__(self):
            ret = ['<real_array>{\n']
            ret.append('    item : ')
            ret.append(repr(self.item))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def testf(self):
        """
        testf(self)
        
        
        Defined at ./test.f90 lines 6-11
        
        Parameters
        ----------
        x : float
        
        """
        _module.f90wrap_testf(x=self._handle)
    
    _dt_array_initialisers = []
    

module_test = Module_Test()

