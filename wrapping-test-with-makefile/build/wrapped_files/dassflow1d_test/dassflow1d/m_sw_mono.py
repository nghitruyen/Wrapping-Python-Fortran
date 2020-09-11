"""
Module m_sw_mono


Defined at m_sw_mono.f90 lines 56-584

"""
from __future__ import print_function, absolute_import, division
import _dassflow1d
import f90wrap.runtime
import logging
from dassflow1d.m_mesh import Mesh

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("dassflow1d.Unknowns")
class Unknowns(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=unknowns)
    
    
    Defined at m_sw_mono.f90 lines 72-84
    
    """
    def __init__(self, msh, handle=None):
        """
        self = Unknowns(msh)
        
        
        Defined at m_sw_mono.f90 lines 554-572
        
        Parameters
        ----------
        msh : Mesh
        
        Returns
        -------
        dof : Unknowns
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_unknowns_initialise(msh=msh._handle)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Unknowns
        
        
        Defined at m_sw_mono.f90 lines 574-583
        
        Parameters
        ----------
        dof : Unknowns
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_unknowns_finalise(dof=self._handle)
    
    @property
    def a(self):
        """
        Element a ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 76
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def qlat(self):
        """
        Element qlat ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 80
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__qlat(self._handle)
        if array_handle in self._arrays:
            qlat = self._arrays[array_handle]
        else:
            qlat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__qlat)
            self._arrays[array_handle] = qlat
        return qlat
    
    @qlat.setter
    def qlat(self, qlat):
        self.qlat[...] = qlat
    
    @property
    def sf(self):
        """
        Element sf ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 82
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__sf(self._handle)
        if array_handle in self._arrays:
            sf = self._arrays[array_handle]
        else:
            sf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__sf)
            self._arrays[array_handle] = sf
        return sf
    
    @sf.setter
    def sf(self, sf):
        self.sf[...] = sf
    
    @property
    def sg(self):
        """
        Element sg ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 84
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_unknowns__array__sg(self._handle)
        if array_handle in self._arrays:
            sg = self._arrays[array_handle]
        else:
            sg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_unknowns__array__sg)
            self._arrays[array_handle] = sg
        return sg
    
    @sg.setter
    def sg(self, sg):
        self.sg[...] = sg
    
    def __str__(self):
        ret = ['<unknowns>{\n']
        ret.append('    a : ')
        ret.append(repr(self.a))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    qlat : ')
        ret.append(repr(self.qlat))
        ret.append(',\n    sf : ')
        ret.append(repr(self.sf))
        ret.append(',\n    sg : ')
        ret.append(repr(self.sg))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.ImplicitMatrix")
class ImplicitMatrix(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=implicitmatrix)
    
    
    Defined at m_sw_mono.f90 lines 87-121
    
    """
    def __init__(self, mdl, msh, handle=None):
        """
        self = Implicitmatrix(mdl, msh)
        
        
        Defined at m_sw_mono.f90 lines 276-312
        
        Parameters
        ----------
        mdl : Model
        msh : Mesh
        
        Returns
        -------
        imp : Implicitmatrix
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_implicitmatrix_initialise(mdl=mdl._handle, \
            msh=msh._handle)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Implicitmatrix
        
        
        Defined at m_sw_mono.f90 lines 314-333
        
        Parameters
        ----------
        imp : Implicitmatrix
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_implicitmatrix_finalise(imp=self._handle)
    
    @property
    def seg_offsets(self):
        """
        Element seg_offsets ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__seg_offsets(self._handle)
        if array_handle in self._arrays:
            seg_offsets = self._arrays[array_handle]
        else:
            seg_offsets = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__seg_offsets)
            self._arrays[array_handle] = seg_offsets
        return seg_offsets
    
    @seg_offsets.setter
    def seg_offsets(self, seg_offsets):
        self.seg_offsets[...] = seg_offsets
    
    @property
    def ga(self):
        """
        Element ga ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 91
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ga(self._handle)
        if array_handle in self._arrays:
            ga = self._arrays[array_handle]
        else:
            ga = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ga)
            self._arrays[array_handle] = ga
        return ga
    
    @ga.setter
    def ga(self, ga):
        self.ga[...] = ga
    
    @property
    def gb(self):
        """
        Element gb ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 93
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__gb(self._handle)
        if array_handle in self._arrays:
            gb = self._arrays[array_handle]
        else:
            gb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__gb)
            self._arrays[array_handle] = gb
        return gb
    
    @gb.setter
    def gb(self, gb):
        self.gb[...] = gb
    
    @property
    def gc(self):
        """
        Element gc ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 95
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__gc(self._handle)
        if array_handle in self._arrays:
            gc = self._arrays[array_handle]
        else:
            gc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__gc)
            self._arrays[array_handle] = gc
        return gc
    
    @gc.setter
    def gc(self, gc):
        self.gc[...] = gc
    
    @property
    def gd(self):
        """
        Element gd ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 97
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__gd(self._handle)
        if array_handle in self._arrays:
            gd = self._arrays[array_handle]
        else:
            gd = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__gd)
            self._arrays[array_handle] = gd
        return gd
    
    @gd.setter
    def gd(self, gd):
        self.gd[...] = gd
    
    @property
    def ge(self):
        """
        Element ge ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 99
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ge(self._handle)
        if array_handle in self._arrays:
            ge = self._arrays[array_handle]
        else:
            ge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ge)
            self._arrays[array_handle] = ge
        return ge
    
    @ge.setter
    def ge(self, ge):
        self.ge[...] = ge
    
    @property
    def gf(self):
        """
        Element gf ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 101
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__gf(self._handle)
        if array_handle in self._arrays:
            gf = self._arrays[array_handle]
        else:
            gf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__gf)
            self._arrays[array_handle] = gf
        return gf
    
    @gf.setter
    def gf(self, gf):
        self.gf[...] = gf
    
    @property
    def cr(self):
        """
        Element cr ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 103
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__cr(self._handle)
        if array_handle in self._arrays:
            cr = self._arrays[array_handle]
        else:
            cr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__cr)
            self._arrays[array_handle] = cr
        return cr
    
    @cr.setter
    def cr(self, cr):
        self.cr[...] = cr
    
    @property
    def cs(self):
        """
        Element cs ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 105
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__cs(self._handle)
        if array_handle in self._arrays:
            cs = self._arrays[array_handle]
        else:
            cs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__cs)
            self._arrays[array_handle] = cs
        return cs
    
    @cs.setter
    def cs(self, cs):
        self.cs[...] = cs
    
    @property
    def ct(self):
        """
        Element ct ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 107
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ct(self._handle)
        if array_handle in self._arrays:
            ct = self._arrays[array_handle]
        else:
            ct = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ct)
            self._arrays[array_handle] = ct
        return ct
    
    @ct.setter
    def ct(self, ct):
        self.ct[...] = ct
    
    @property
    def ta1(self):
        """
        Element ta1 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 109
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ta1(self._handle)
        if array_handle in self._arrays:
            ta1 = self._arrays[array_handle]
        else:
            ta1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ta1)
            self._arrays[array_handle] = ta1
        return ta1
    
    @ta1.setter
    def ta1(self, ta1):
        self.ta1[...] = ta1
    
    @property
    def ta2(self):
        """
        Element ta2 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 110
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ta2(self._handle)
        if array_handle in self._arrays:
            ta2 = self._arrays[array_handle]
        else:
            ta2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ta2)
            self._arrays[array_handle] = ta2
        return ta2
    
    @ta2.setter
    def ta2(self, ta2):
        self.ta2[...] = ta2
    
    @property
    def ta3(self):
        """
        Element ta3 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 111
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ta3(self._handle)
        if array_handle in self._arrays:
            ta3 = self._arrays[array_handle]
        else:
            ta3 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ta3)
            self._arrays[array_handle] = ta3
        return ta3
    
    @ta3.setter
    def ta3(self, ta3):
        self.ta3[...] = ta3
    
    @property
    def ta4(self):
        """
        Element ta4 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 112
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__ta4(self._handle)
        if array_handle in self._arrays:
            ta4 = self._arrays[array_handle]
        else:
            ta4 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__ta4)
            self._arrays[array_handle] = ta4
        return ta4
    
    @ta4.setter
    def ta4(self, ta4):
        self.ta4[...] = ta4
    
    @property
    def tb1(self):
        """
        Element tb1 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 114
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__tb1(self._handle)
        if array_handle in self._arrays:
            tb1 = self._arrays[array_handle]
        else:
            tb1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__tb1)
            self._arrays[array_handle] = tb1
        return tb1
    
    @tb1.setter
    def tb1(self, tb1):
        self.tb1[...] = tb1
    
    @property
    def tb2(self):
        """
        Element tb2 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 115
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__tb2(self._handle)
        if array_handle in self._arrays:
            tb2 = self._arrays[array_handle]
        else:
            tb2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__tb2)
            self._arrays[array_handle] = tb2
        return tb2
    
    @tb2.setter
    def tb2(self, tb2):
        self.tb2[...] = tb2
    
    @property
    def anz(self):
        """
        Element anz ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 120
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__anz(self._handle)
        if array_handle in self._arrays:
            anz = self._arrays[array_handle]
        else:
            anz = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__anz)
            self._arrays[array_handle] = anz
        return anz
    
    @anz.setter
    def anz(self, anz):
        self.anz[...] = anz
    
    @property
    def rhs(self):
        """
        Element rhs ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 121
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_implicitmatrix__array__rhs(self._handle)
        if array_handle in self._arrays:
            rhs = self._arrays[array_handle]
        else:
            rhs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_implicitmatrix__array__rhs)
            self._arrays[array_handle] = rhs
        return rhs
    
    @rhs.setter
    def rhs(self, rhs):
        self.rhs[...] = rhs
    
    def __str__(self):
        ret = ['<implicitmatrix>{\n']
        ret.append('    seg_offsets : ')
        ret.append(repr(self.seg_offsets))
        ret.append(',\n    ga : ')
        ret.append(repr(self.ga))
        ret.append(',\n    gb : ')
        ret.append(repr(self.gb))
        ret.append(',\n    gc : ')
        ret.append(repr(self.gc))
        ret.append(',\n    gd : ')
        ret.append(repr(self.gd))
        ret.append(',\n    ge : ')
        ret.append(repr(self.ge))
        ret.append(',\n    gf : ')
        ret.append(repr(self.gf))
        ret.append(',\n    cr : ')
        ret.append(repr(self.cr))
        ret.append(',\n    cs : ')
        ret.append(repr(self.cs))
        ret.append(',\n    ct : ')
        ret.append(repr(self.ct))
        ret.append(',\n    ta1 : ')
        ret.append(repr(self.ta1))
        ret.append(',\n    ta2 : ')
        ret.append(repr(self.ta2))
        ret.append(',\n    ta3 : ')
        ret.append(repr(self.ta3))
        ret.append(',\n    ta4 : ')
        ret.append(repr(self.ta4))
        ret.append(',\n    tb1 : ')
        ret.append(repr(self.tb1))
        ret.append(',\n    tb2 : ')
        ret.append(repr(self.tb2))
        ret.append(',\n    anz : ')
        ret.append(repr(self.anz))
        ret.append(',\n    rhs : ')
        ret.append(repr(self.rhs))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Timeseries")
class Timeseries(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=timeseries)
    
    
    Defined at m_sw_mono.f90 lines 124-126
    
    """
    def __init__(self, nt, y0=None, handle=None):
        """
        self = Timeseries(nt[, y0])
        
        
        Defined at m_sw_mono.f90 lines 529-544
        
        Parameters
        ----------
        nt : int
        y0 : float
        
        Returns
        -------
        ts : Timeseries
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_timeseries_initialise(nt=nt, y0=y0)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Timeseries
        
        
        Defined at m_sw_mono.f90 lines 546-551
        
        Parameters
        ----------
        ts : Timeseries
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_timeseries_finalise(ts=self._handle)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 125
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_timeseries__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_timeseries__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 126
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_timeseries__array__y(self._handle)
        if array_handle in self._arrays:
            y = self._arrays[array_handle]
        else:
            y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_timeseries__array__y)
            self._arrays[array_handle] = y
        return y
    
    @y.setter
    def y(self, y):
        self.y[...] = y
    
    def __str__(self):
        ret = ['<timeseries>{\n']
        ret.append('    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.BoundaryCondition")
class BoundaryCondition(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=boundarycondition)
    
    
    Defined at m_sw_mono.f90 lines 129-133
    
    """
    def set_timeseries(self, t, y):
        """
        set_timeseries(self, t, y)
        
        
        Defined at m_sw_mono.f90 lines 261-274
        
        Parameters
        ----------
        bc : Boundarycondition
        t : float array
        y : float array
        
        """
        _dassflow1d.f90wrap_set_timeseries(bc=self._handle, t=t, y=y)
    
    def __init__(self, handle=None):
        """
        self = Boundarycondition()
        
        
        Defined at m_sw_mono.f90 lines 129-133
        
        
        Returns
        -------
        this : Boundarycondition
        	Object to be constructed
        
        
        Automatically generated constructor for boundarycondition
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_boundarycondition_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Boundarycondition
        
        
        Defined at m_sw_mono.f90 lines 129-133
        
        Parameters
        ----------
        this : Boundarycondition
        	Object to be destructed
        
        
        Automatically generated destructor for boundarycondition
        """
        if self._alloc:
            _dassflow1d.f90wrap_boundarycondition_finalise(this=self._handle)
    
    @property
    def id(self):
        """
        Element id ftype=character(len=16) pytype=str
        
        
        Defined at m_sw_mono.f90 line 131
        
        """
        return _dassflow1d.f90wrap_boundarycondition__get__id(self._handle)
    
    @id.setter
    def id(self, id):
        _dassflow1d.f90wrap_boundarycondition__set__id(self._handle, id)
    
    @property
    def ts(self):
        """
        Element ts ftype=type(timeseries) pytype=Timeseries
        
        
        Defined at m_sw_mono.f90 line 133
        
        """
        ts_handle = _dassflow1d.f90wrap_boundarycondition__get__ts(self._handle)
        if tuple(ts_handle) in self._objs:
            ts = self._objs[tuple(ts_handle)]
        else:
            ts = Timeseries.from_handle(ts_handle)
            self._objs[tuple(ts_handle)] = ts
        return ts
    
    @ts.setter
    def ts(self, ts):
        ts = ts._handle
        _dassflow1d.f90wrap_boundarycondition__set__ts(self._handle, ts)
    
    def __str__(self):
        ret = ['<boundarycondition>{\n']
        ret.append('    id : ')
        ret.append(repr(self.id))
        ret.append(',\n    ts : ')
        ret.append(repr(self.ts))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.InflowCondition")
class InflowCondition(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=inflowcondition)
    
    
    Defined at m_sw_mono.f90 lines 136-142
    
    """
    def __init__(self, handle=None):
        """
        self = Inflowcondition()
        
        
        Defined at m_sw_mono.f90 lines 136-142
        
        
        Returns
        -------
        this : Inflowcondition
        	Object to be constructed
        
        
        Automatically generated constructor for inflowcondition
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_inflowcondition_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Inflowcondition
        
        
        Defined at m_sw_mono.f90 lines 136-142
        
        Parameters
        ----------
        this : Inflowcondition
        	Object to be destructed
        
        
        Automatically generated destructor for inflowcondition
        """
        if self._alloc:
            _dassflow1d.f90wrap_inflowcondition_finalise(this=self._handle)
    
    @property
    def iseg(self):
        """
        Element iseg ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 138
        
        """
        return _dassflow1d.f90wrap_inflowcondition__get__iseg(self._handle)
    
    @iseg.setter
    def iseg(self, iseg):
        _dassflow1d.f90wrap_inflowcondition__set__iseg(self._handle, iseg)
    
    @property
    def ie(self):
        """
        Element ie ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 140
        
        """
        return _dassflow1d.f90wrap_inflowcondition__get__ie(self._handle)
    
    @ie.setter
    def ie(self, ie):
        _dassflow1d.f90wrap_inflowcondition__set__ie(self._handle, ie)
    
    @property
    def ts(self):
        """
        Element ts ftype=type(timeseries) pytype=Timeseries
        
        
        Defined at m_sw_mono.f90 line 142
        
        """
        ts_handle = _dassflow1d.f90wrap_inflowcondition__get__ts(self._handle)
        if tuple(ts_handle) in self._objs:
            ts = self._objs[tuple(ts_handle)]
        else:
            ts = Timeseries.from_handle(ts_handle)
            self._objs[tuple(ts_handle)] = ts
        return ts
    
    @ts.setter
    def ts(self, ts):
        ts = ts._handle
        _dassflow1d.f90wrap_inflowcondition__set__ts(self._handle, ts)
    
    def __str__(self):
        ret = ['<inflowcondition>{\n']
        ret.append('    iseg : ')
        ret.append(repr(self.iseg))
        ret.append(',\n    ie : ')
        ret.append(repr(self.ie))
        ret.append(',\n    ts : ')
        ret.append(repr(self.ts))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Results")
class Results(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=results)
    
    
    Defined at m_sw_mono.f90 lines 144-148
    
    """
    def __init__(self, handle=None):
        """
        self = Results()
        
        
        Defined at m_sw_mono.f90 lines 144-148
        
        
        Returns
        -------
        this : Results
        	Object to be constructed
        
        
        Automatically generated constructor for results
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_results_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Results
        
        
        Defined at m_sw_mono.f90 lines 144-148
        
        Parameters
        ----------
        this : Results
        	Object to be destructed
        
        
        Automatically generated destructor for results
        """
        if self._alloc:
            _dassflow1d.f90wrap_results_finalise(this=self._handle)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 145
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_results__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_results__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 146
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_results__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_results__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 147
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_results__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_results__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def a(self):
        """
        Element a ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_results__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_results__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    def __str__(self):
        ret = ['<results>{\n']
        ret.append('    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    a : ')
        ret.append(repr(self.a))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Model")
class Model(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=model)
    
    
    Defined at m_sw_mono.f90 lines 150-203
    
    """
    def __init__(self, msh, scheme=None, handle=None):
        """
        self = Model(msh[, scheme])
        
        
        Defined at m_sw_mono.f90 lines 335-383
        
        Parameters
        ----------
        msh : Mesh
        scheme : str
        
        Returns
        -------
        mdl : Model
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_model_initialise(msh=msh._handle, scheme=scheme)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def set_large_grid(self, large_grid):
        """
        set_large_grid(self, large_grid)
        
        
        Defined at m_sw_mono.f90 lines 385-389
        
        Parameters
        ----------
        mdl : Model
        large_grid : Mesh
        
        """
        _dassflow1d.f90wrap_set_large_grid(mdl=self._handle, \
            large_grid=large_grid._handle)
    
    def add_inflow_condition(self, iseg, coords, t, q):
        """
        add_inflow_condition(self, iseg, coords, t, q)
        
        
        Defined at m_sw_mono.f90 lines 391-486
        
        Parameters
        ----------
        mdl : Model
        iseg : int
        coords : float array
        t : float array
        q : float array
        
        """
        _dassflow1d.f90wrap_add_inflow_condition(mdl=self._handle, iseg=iseg, \
            coords=coords, t=t, q=q)
    
    def set_scheme(self, scheme=None):
        """
        set_scheme(self[, scheme])
        
        
        Defined at m_sw_mono.f90 lines 488-518
        
        Parameters
        ----------
        mdl : Model
        scheme : str
        
        """
        _dassflow1d.f90wrap_set_scheme(mdl=self._handle, scheme=scheme)
    
    def __del__(self):
        """
        Destructor for class Model
        
        
        Defined at m_sw_mono.f90 lines 520-527
        
        Parameters
        ----------
        mdl : Model
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_model_finalise(mdl=self._handle)
    
    @property
    def discharge_estimation(self):
        """
        Element discharge_estimation ftype=logical pytype=bool
        
        
        Defined at m_sw_mono.f90 line 152
        
        """
        return _dassflow1d.f90wrap_model__get__discharge_estimation(self._handle)
    
    @discharge_estimation.setter
    def discharge_estimation(self, discharge_estimation):
        _dassflow1d.f90wrap_model__set__discharge_estimation(self._handle, \
            discharge_estimation)
    
    @property
    def status(self):
        """
        Element status ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 154
        
        """
        return _dassflow1d.f90wrap_model__get__status(self._handle)
    
    @status.setter
    def status(self, status):
        _dassflow1d.f90wrap_model__set__status(self._handle, status)
    
    @property
    def warning_counters(self):
        """
        Element warning_counters ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 156
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_model__array__warning_counters(self._handle)
        if array_handle in self._arrays:
            warning_counters = self._arrays[array_handle]
        else:
            warning_counters = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_model__array__warning_counters)
            self._arrays[array_handle] = warning_counters
        return warning_counters
    
    @warning_counters.setter
    def warning_counters(self, warning_counters):
        self.warning_counters[...] = warning_counters
    
    @property
    def tc(self):
        """
        Element tc ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 158
        
        """
        return _dassflow1d.f90wrap_model__get__tc(self._handle)
    
    @tc.setter
    def tc(self, tc):
        _dassflow1d.f90wrap_model__set__tc(self._handle, tc)
    
    @property
    def te(self):
        """
        Element te ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 160
        
        """
        return _dassflow1d.f90wrap_model__get__te(self._handle)
    
    @te.setter
    def te(self, te):
        _dassflow1d.f90wrap_model__set__te(self._handle, te)
    
    @property
    def ts(self):
        """
        Element ts ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 162
        
        """
        return _dassflow1d.f90wrap_model__get__ts(self._handle)
    
    @ts.setter
    def ts(self, ts):
        _dassflow1d.f90wrap_model__set__ts(self._handle, ts)
    
    @property
    def dt(self):
        """
        Element dt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 164
        
        """
        return _dassflow1d.f90wrap_model__get__dt(self._handle)
    
    @dt.setter
    def dt(self, dt):
        _dassflow1d.f90wrap_model__set__dt(self._handle, dt)
    
    @property
    def dtout(self):
        """
        Element dtout ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 166
        
        """
        return _dassflow1d.f90wrap_model__get__dtout(self._handle)
    
    @dtout.setter
    def dtout(self, dtout):
        _dassflow1d.f90wrap_model__set__dtout(self._handle, dtout)
    
    @property
    def frlpi(self):
        """
        Element frlpi ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 168
        
        """
        return _dassflow1d.f90wrap_model__get__frlpi(self._handle)
    
    @frlpi.setter
    def frlpi(self, frlpi):
        _dassflow1d.f90wrap_model__set__frlpi(self._handle, frlpi)
    
    @property
    def gravity(self):
        """
        Element gravity ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 170
        
        """
        return _dassflow1d.f90wrap_model__get__gravity(self._handle)
    
    @gravity.setter
    def gravity(self, gravity):
        _dassflow1d.f90wrap_model__set__gravity(self._handle, gravity)
    
    @property
    def heps(self):
        """
        Element heps ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 172
        
        """
        return _dassflow1d.f90wrap_model__get__heps(self._handle)
    
    @heps.setter
    def heps(self, heps):
        _dassflow1d.f90wrap_model__set__heps(self._handle, heps)
    
    @property
    def qeps(self):
        """
        Element qeps ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 174
        
        """
        return _dassflow1d.f90wrap_model__get__qeps(self._handle)
    
    @qeps.setter
    def qeps(self, qeps):
        _dassflow1d.f90wrap_model__set__qeps(self._handle, qeps)
    
    @property
    def mlpi(self):
        """
        Element mlpi ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 176
        
        """
        return _dassflow1d.f90wrap_model__get__mlpi(self._handle)
    
    @mlpi.setter
    def mlpi(self, mlpi):
        _dassflow1d.f90wrap_model__set__mlpi(self._handle, mlpi)
    
    @property
    def theta_preissmann(self):
        """
        Element theta_preissmann ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 178
        
        """
        return _dassflow1d.f90wrap_model__get__theta_preissmann(self._handle)
    
    @theta_preissmann.setter
    def theta_preissmann(self, theta_preissmann):
        _dassflow1d.f90wrap_model__set__theta_preissmann(self._handle, theta_preissmann)
    
    @property
    def gamma_reg(self):
        """
        Element gamma_reg ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 181
        
        """
        return _dassflow1d.f90wrap_model__get__gamma_reg(self._handle)
    
    @gamma_reg.setter
    def gamma_reg(self, gamma_reg):
        _dassflow1d.f90wrap_model__set__gamma_reg(self._handle, gamma_reg)
    
    @property
    def cost_obs(self):
        """
        Element cost_obs ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        return _dassflow1d.f90wrap_model__get__cost_obs(self._handle)
    
    @cost_obs.setter
    def cost_obs(self, cost_obs):
        _dassflow1d.f90wrap_model__set__cost_obs(self._handle, cost_obs)
    
    @property
    def cost_reg(self):
        """
        Element cost_reg ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 185
        
        """
        return _dassflow1d.f90wrap_model__get__cost_reg(self._handle)
    
    @cost_reg.setter
    def cost_reg(self, cost_reg):
        _dassflow1d.f90wrap_model__set__cost_reg(self._handle, cost_reg)
    
    @property
    def scheme(self):
        """
        Element scheme ftype=character(len=32) pytype=str
        
        
        Defined at m_sw_mono.f90 line 187
        
        """
        return _dassflow1d.f90wrap_model__get__scheme(self._handle)
    
    @scheme.setter
    def scheme(self, scheme):
        _dassflow1d.f90wrap_model__set__scheme(self._handle, scheme)
    
    @property
    def scale_idw(self):
        """
        Element scale_idw ftype=character(len=16) pytype=str
        
        
        Defined at m_sw_mono.f90 line 189
        
        """
        return _dassflow1d.f90wrap_model__get__scale_idw(self._handle)
    
    @scale_idw.setter
    def scale_idw(self, scale_idw):
        _dassflow1d.f90wrap_model__set__scale_idw(self._handle, scale_idw)
    
    def init_array_bc(self):
        self.bc = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_model__array_getitem__bc,
                                        _dassflow1d.f90wrap_model__array_setitem__bc,
                                        _dassflow1d.f90wrap_model__array_len__bc,
                                        """
        Element bc ftype=type(boundarycondition) pytype=Boundarycondition
        
        
        Defined at m_sw_mono.f90 line 191
        
        """, BoundaryCondition)
        return self.bc
    
    def init_array_ic(self):
        self.ic = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_model__array_getitem__ic,
                                        _dassflow1d.f90wrap_model__array_setitem__ic,
                                        _dassflow1d.f90wrap_model__array_len__ic,
                                        """
        Element ic ftype=type(inflowcondition) pytype=Inflowcondition
        
        
        Defined at m_sw_mono.f90 line 193
        
        """, InflowCondition)
        return self.ic
    
    @property
    def msh(self):
        """
        Element msh ftype=type(mesh) pytype=Mesh
        
        
        Defined at m_sw_mono.f90 line 195
        
        """
        msh_handle = _dassflow1d.f90wrap_model__get__msh(self._handle)
        if tuple(msh_handle) in self._objs:
            msh = self._objs[tuple(msh_handle)]
        else:
            msh = Mesh.from_handle(msh_handle)
            self._objs[tuple(msh_handle)] = msh
        return msh
    
    @msh.setter
    def msh(self, msh):
        msh = msh._handle
        _dassflow1d.f90wrap_model__set__msh(self._handle, msh)
    
    @property
    def large_grid(self):
        """
        Element large_grid ftype=type(mesh) pytype=Mesh
        
        
        Defined at m_sw_mono.f90 line 197
        
        """
        large_grid_handle = _dassflow1d.f90wrap_model__get__large_grid(self._handle)
        if tuple(large_grid_handle) in self._objs:
            large_grid = self._objs[tuple(large_grid_handle)]
        else:
            large_grid = Mesh.from_handle(large_grid_handle)
            self._objs[tuple(large_grid_handle)] = large_grid
        return large_grid
    
    @large_grid.setter
    def large_grid(self, large_grid):
        large_grid = large_grid._handle
        _dassflow1d.f90wrap_model__set__large_grid(self._handle, large_grid)
    
    @property
    def dof(self):
        """
        Element dof ftype=type(unknowns) pytype=Unknowns
        
        
        Defined at m_sw_mono.f90 line 199
        
        """
        dof_handle = _dassflow1d.f90wrap_model__get__dof(self._handle)
        if tuple(dof_handle) in self._objs:
            dof = self._objs[tuple(dof_handle)]
        else:
            dof = Unknowns.from_handle(dof_handle)
            self._objs[tuple(dof_handle)] = dof
        return dof
    
    @dof.setter
    def dof(self, dof):
        dof = dof._handle
        _dassflow1d.f90wrap_model__set__dof(self._handle, dof)
    
    @property
    def imp(self):
        """
        Element imp ftype=type(implicitmatrix) pytype=Implicitmatrix
        
        
        Defined at m_sw_mono.f90 line 201
        
        """
        imp_handle = _dassflow1d.f90wrap_model__get__imp(self._handle)
        if tuple(imp_handle) in self._objs:
            imp = self._objs[tuple(imp_handle)]
        else:
            imp = ImplicitMatrix.from_handle(imp_handle)
            self._objs[tuple(imp_handle)] = imp
        return imp
    
    @imp.setter
    def imp(self, imp):
        imp = imp._handle
        _dassflow1d.f90wrap_model__set__imp(self._handle, imp)
    
    @property
    def res(self):
        """
        Element res ftype=type(results) pytype=Results
        
        
        Defined at m_sw_mono.f90 line 203
        
        """
        res_handle = _dassflow1d.f90wrap_model__get__res(self._handle)
        if tuple(res_handle) in self._objs:
            res = self._objs[tuple(res_handle)]
        else:
            res = Results.from_handle(res_handle)
            self._objs[tuple(res_handle)] = res
        return res
    
    @res.setter
    def res(self, res):
        res = res._handle
        _dassflow1d.f90wrap_model__set__res(self._handle, res)
    
    def __str__(self):
        ret = ['<model>{\n']
        ret.append('    discharge_estimation : ')
        ret.append(repr(self.discharge_estimation))
        ret.append(',\n    status : ')
        ret.append(repr(self.status))
        ret.append(',\n    warning_counters : ')
        ret.append(repr(self.warning_counters))
        ret.append(',\n    tc : ')
        ret.append(repr(self.tc))
        ret.append(',\n    te : ')
        ret.append(repr(self.te))
        ret.append(',\n    ts : ')
        ret.append(repr(self.ts))
        ret.append(',\n    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    dtout : ')
        ret.append(repr(self.dtout))
        ret.append(',\n    frlpi : ')
        ret.append(repr(self.frlpi))
        ret.append(',\n    gravity : ')
        ret.append(repr(self.gravity))
        ret.append(',\n    heps : ')
        ret.append(repr(self.heps))
        ret.append(',\n    qeps : ')
        ret.append(repr(self.qeps))
        ret.append(',\n    mlpi : ')
        ret.append(repr(self.mlpi))
        ret.append(',\n    theta_preissmann : ')
        ret.append(repr(self.theta_preissmann))
        ret.append(',\n    gamma_reg : ')
        ret.append(repr(self.gamma_reg))
        ret.append(',\n    cost_obs : ')
        ret.append(repr(self.cost_obs))
        ret.append(',\n    cost_reg : ')
        ret.append(repr(self.cost_reg))
        ret.append(',\n    scheme : ')
        ret.append(repr(self.scheme))
        ret.append(',\n    scale_idw : ')
        ret.append(repr(self.scale_idw))
        ret.append(',\n    msh : ')
        ret.append(repr(self.msh))
        ret.append(',\n    large_grid : ')
        ret.append(repr(self.large_grid))
        ret.append(',\n    dof : ')
        ret.append(repr(self.dof))
        ret.append(',\n    imp : ')
        ret.append(repr(self.imp))
        ret.append(',\n    res : ')
        ret.append(repr(self.res))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_bc, init_array_ic]
    


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_sw_mono".')

for func in _dt_array_initialisers:
    func()
