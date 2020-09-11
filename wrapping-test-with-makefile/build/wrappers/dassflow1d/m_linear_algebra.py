"""
Module m_linear_algebra


Defined at m_linear_algebra.f90 lines 57-650

"""
import _dassflow1d
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

class Vec2D(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=vec2d)
    
    
    Defined at m_linear_algebra.f90 lines 64-69
    
    """
    def __init__(self, handle=None):
        """
        self = Vec2D()
        
        
        Defined at m_linear_algebra.f90 lines 64-69
        
        
        Returns
        -------
        this : Vec2D
        	Object to be constructed
        
        
        Automatically generated constructor for vec2d
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _dassflow1d.f90wrap_vec2d_initialise()
    
    def __del__(self):
        """
        Destructor for class Vec2D
        
        
        Defined at m_linear_algebra.f90 lines 64-69
        
        Parameters
        ----------
        this : Vec2D
        	Object to be destructed
        
        
        Automatically generated destructor for vec2d
        """
        if self._alloc:
            _dassflow1d.f90wrap_vec2d_finalise(this=self._handle)
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 66
        
        """
        return _dassflow1d.f90wrap_vec2d__get__x(self._handle)
    
    @x.setter
    def x(self, x):
        _dassflow1d.f90wrap_vec2d__set__x(self._handle, x)
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 68
        
        """
        return _dassflow1d.f90wrap_vec2d__get__y(self._handle)
    
    @y.setter
    def y(self, y):
        _dassflow1d.f90wrap_vec2d__set__y(self._handle, y)
    
    def __str__(self):
        ret = ['<vec2d>{\n']
        ret.append('    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    
class Matrix(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=matrix)
    
    
    Defined at m_linear_algebra.f90 lines 72-76
    
    """
    def __init__(self, handle=None):
        """
        self = Matrix()
        
        
        Defined at m_linear_algebra.f90 lines 72-76
        
        
        Returns
        -------
        this : Matrix
        	Object to be constructed
        
        
        Automatically generated constructor for matrix
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _dassflow1d.f90wrap_matrix_initialise()
    
    def __del__(self):
        """
        Destructor for class Matrix
        
        
        Defined at m_linear_algebra.f90 lines 72-76
        
        Parameters
        ----------
        this : Matrix
        	Object to be destructed
        
        
        Automatically generated destructor for matrix
        """
        if self._alloc:
            _dassflow1d.f90wrap_matrix_finalise(this=self._handle)
    
    @property
    def m(self):
        """
        Element m ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_matrix__array__m(self._handle)
        if array_handle in self._arrays:
            m = self._arrays[array_handle]
        else:
            m = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_matrix__array__m)
            self._arrays[array_handle] = m
        return m
    
    @m.setter
    def m(self, m):
        self.m[...] = m
    
    def __str__(self):
        ret = ['<matrix>{\n']
        ret.append('    m : ')
        ret.append(repr(self.m))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    
class Matrixcsr(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=matrixcsr)
    
    
    Defined at m_linear_algebra.f90 lines 79-91
    
    """
    def __init__(self, n, nnz, handle=None):
        """
        self = Matrixcsr(n, nnz)
        
        
        Defined at m_linear_algebra.f90 lines 139-153
        
        Parameters
        ----------
        n : int
        nnz : int
        
        Returns
        -------
        mat : Matrixcsr
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _dassflow1d.f90wrap_matrixcsr_initialise(n=n, nnz=nnz)
    
    def csr_x_vec_add(self, x, y):
        """
        csr_x_vec_add(self, x, y)
        
        
        Defined at m_linear_algebra.f90 lines 626-650
        
        Parameters
        ----------
        a : Matrixcsr
        x : float array
        y : float array
        
        """
        _dassflow1d.f90wrap_csr_x_vec_add(a=self._handle, x=x, y=y)
    
    def __del__(self):
        """
        Destructor for class Matrixcsr
        
        
        Defined at m_linear_algebra.f90 lines 79-91
        
        Parameters
        ----------
        this : Matrixcsr
        	Object to be destructed
        
        
        Automatically generated destructor for matrixcsr
        """
        if self._alloc:
            _dassflow1d.f90wrap_matrixcsr_finalise(this=self._handle)
    
    @property
    def n(self):
        """
        Element n ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 81
        
        """
        return _dassflow1d.f90wrap_matrixcsr__get__n(self._handle)
    
    @n.setter
    def n(self, n):
        _dassflow1d.f90wrap_matrixcsr__set__n(self._handle, n)
    
    @property
    def nnz(self):
        """
        Element nnz ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 83
        
        """
        return _dassflow1d.f90wrap_matrixcsr__get__nnz(self._handle)
    
    @nnz.setter
    def nnz(self, nnz):
        _dassflow1d.f90wrap_matrixcsr__set__nnz(self._handle, nnz)
    
    @property
    def irow(self):
        """
        Element irow ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 85
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_matrixcsr__array__irow(self._handle)
        if array_handle in self._arrays:
            irow = self._arrays[array_handle]
        else:
            irow = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_matrixcsr__array__irow)
            self._arrays[array_handle] = irow
        return irow
    
    @irow.setter
    def irow(self, irow):
        self.irow[...] = irow
    
    @property
    def icol(self):
        """
        Element icol ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 87
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_matrixcsr__array__icol(self._handle)
        if array_handle in self._arrays:
            icol = self._arrays[array_handle]
        else:
            icol = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_matrixcsr__array__icol)
            self._arrays[array_handle] = icol
        return icol
    
    @icol.setter
    def icol(self, icol):
        self.icol[...] = icol
    
    @property
    def anz(self):
        """
        Element anz ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_matrixcsr__array__anz(self._handle)
        if array_handle in self._arrays:
            anz = self._arrays[array_handle]
        else:
            anz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_matrixcsr__array__anz)
            self._arrays[array_handle] = anz
        return anz
    
    @anz.setter
    def anz(self, anz):
        self.anz[...] = anz
    
    def __str__(self):
        ret = ['<matrixcsr>{\n']
        ret.append('    n : ')
        ret.append(repr(self.n))
        ret.append(',\n    nnz : ')
        ret.append(repr(self.nnz))
        ret.append(',\n    irow : ')
        ret.append(repr(self.irow))
        ret.append(',\n    icol : ')
        ret.append(repr(self.icol))
        ret.append(',\n    anz : ')
        ret.append(repr(self.anz))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    
class Matrixblock(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=matrixblock)
    
    
    Defined at m_linear_algebra.f90 lines 94-101
    
    """
    def __init__(self, n, sizes, handle=None):
        """
        self = Matrixblock(n, sizes)
        
        
        Defined at m_linear_algebra.f90 lines 107-121
        
        Parameters
        ----------
        n : int
        sizes : int array
        
        Returns
        -------
        mat : Matrixblock
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        self._handle = _dassflow1d.f90wrap_matrixblock_initialise(n=n, sizes=sizes)
    
    def matrixblock_copy(self):
        """
        dst = matrixblock_copy(self)
        
        
        Defined at m_linear_algebra.f90 lines 123-137
        
        Parameters
        ----------
        src : Matrixblock
        
        Returns
        -------
        dst : Matrixblock
        
        """
        dst = _dassflow1d.f90wrap_matrixblock_copy(src=self._handle)
        dst = Matrixblock.from_handle(dst)
        return dst
    
    def matrixcsr_from_matrixblock(self, matcsr, threshold=None):
        """
        matrixcsr_from_matrixblock(self, matcsr[, threshold])
        
        
        Defined at m_linear_algebra.f90 lines 155-218
        
        Parameters
        ----------
        matblock : Matrixblock
        matcsr : Matrixcsr
        threshold : float
        
        """
        _dassflow1d.f90wrap_matrixcsr_from_matrixblock(matblock=self._handle, \
            matcsr=matcsr._handle, threshold=threshold)
    
    def __del__(self):
        """
        Destructor for class Matrixblock
        
        
        Defined at m_linear_algebra.f90 lines 94-101
        
        Parameters
        ----------
        this : Matrixblock
        	Object to be destructed
        
        
        Automatically generated destructor for matrixblock
        """
        if self._alloc:
            _dassflow1d.f90wrap_matrixblock_finalise(this=self._handle)
    
    @property
    def n(self):
        """
        Element n ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 96
        
        """
        return _dassflow1d.f90wrap_matrixblock__get__n(self._handle)
    
    @n.setter
    def n(self, n):
        _dassflow1d.f90wrap_matrixblock__set__n(self._handle, n)
    
    def init_array_blocks(self):
        self.blocks = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_matrixblock__array_getitem__blocks,
                                        _dassflow1d.f90wrap_matrixblock__array_setitem__blocks,
                                        _dassflow1d.f90wrap_matrixblock__array_len__blocks,
                                        """
        Element blocks ftype=type(matrix) pytype=Matrix
        
        
        Defined at m_linear_algebra.f90 line 98
        
        """, Matrix)
        return self.blocks
    
    def __str__(self):
        ret = ['<matrixblock>{\n']
        ret.append('    n : ')
        ret.append(repr(self.n))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_blocks]
    
def matrixcsr_from_numpy_array(array, threshold=None):
    """
    mat = matrixcsr_from_numpy_array(array[, threshold])
    
    
    Defined at m_linear_algebra.f90 lines 220-268
    
    Parameters
    ----------
    array : float array
    threshold : float
    
    Returns
    -------
    mat : Matrixcsr
    
    """
    mat = _dassflow1d.f90wrap_matrixcsr_from_numpy_array(array=array, \
        threshold=threshold)
    mat = Matrixcsr.from_handle(mat)
    return mat

def cholesky_inplace(n, a):
    """
    cholesky_inplace(n, a)
    
    
    Defined at m_linear_algebra.f90 lines 373-467
    
    Parameters
    ----------
    n : int
    a : float array
    
    \
        ================================================================================================================
      Interface Variables
    \
        ================================================================================================================
    """
    _dassflow1d.f90wrap_cholesky_inplace(n=n, a=a)

def solve_using_cholesky_inplace(n, mat, b, x):
    """
    solve_using_cholesky_inplace(n, mat, b, x)
    
    
    Defined at m_linear_algebra.f90 lines 566-621
    
    Parameters
    ----------
    n : int
    mat : Matrix
    b : float array
    x : float array
    
    \
        ================================================================================================================
      Interface Variables
    \
        ================================================================================================================
    """
    _dassflow1d.f90wrap_solve_using_cholesky_inplace(n=n, mat=mat._handle, b=b, x=x)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "m_linear_algebra".')

for func in _dt_array_initialisers:
    func()
