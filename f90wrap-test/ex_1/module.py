from __future__ import print_function, absolute_import, division
import _module
import f90wrap.runtime
import logging

def foo(a, b):
    """
    foo(a, b)
    
    
    Defined at test.f90 lines 1-3
    
    Parameters
    ----------
    a : float
    b : int
    
    """
    _module.f90wrap_foo(a=a, b=b)

