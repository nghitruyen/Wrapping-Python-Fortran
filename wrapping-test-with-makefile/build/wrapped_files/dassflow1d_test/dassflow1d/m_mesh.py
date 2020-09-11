"""
Module m_mesh


Defined at m_mesh.f90 lines 68-827

"""
from __future__ import print_function, absolute_import, division
import _dassflow1d
import f90wrap.runtime
import logging
from dassflow1d.m_linear_algebra import vec2d

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("dassflow1d.Crosssection")
class Crosssection(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=crosssection)
    
    
    Defined at m_mesh.f90 lines 77-114
    
    """
    def __init__(self, nlevels, shape_model=None, handle=None):
        """
        self = Crosssection(nlevels[, shape_model])
        
        
        Defined at m_mesh.f90 lines 373-412
        
        Parameters
        ----------
        nlevels : int
        shape_model : str
        
        Returns
        -------
        cs : Crosssection
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_crosssection_initialise(nlevels=nlevels, \
            shape_model=shape_model)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Crosssection
        
        
        Defined at m_mesh.f90 lines 414-422
        
        Parameters
        ----------
        cs : Crosssection
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_crosssection_finalise(cs=self._handle)
    
    def set_levels(self, heights, widths, shape_model=None):
        """
        set_levels(self, heights, widths[, shape_model])
        
        
        Defined at m_mesh.f90 lines 424-462
        
        Parameters
        ----------
        cs : Crosssection
        heights : float array
        widths : float array
        shape_model : str
        
        """
        _dassflow1d.f90wrap_set_levels(cs=self._handle, heights=heights, widths=widths, \
            shape_model=shape_model)
    
    def crosssection_copy(self, cs_dst):
        """
        crosssection_copy(self, cs_dst)
        
        
        Defined at m_mesh.f90 lines 464-511
        
        Parameters
        ----------
        cs_src : Crosssection
        cs_dst : Crosssection
        
        """
        _dassflow1d.f90wrap_crosssection_copy(cs_src=self._handle, \
            cs_dst=cs_dst._handle)
    
    def update_geometry(self):
        """
        update_geometry(self)
        
        
        Defined at m_mesh.f90 lines 513-517
        
        Parameters
        ----------
        cs : Crosssection
        
        """
        _dassflow1d.f90wrap_update_geometry(cs=self._handle)
    
    @property
    def coord(self):
        """
        Element coord ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 79
        
        """
        coord_handle = _dassflow1d.f90wrap_crosssection__get__coord(self._handle)
        if tuple(coord_handle) in self._objs:
            coord = self._objs[tuple(coord_handle)]
        else:
            coord = vec2d.from_handle(coord_handle)
            self._objs[tuple(coord_handle)] = coord
        return coord
    
    @coord.setter
    def coord(self, coord):
        coord = coord._handle
        _dassflow1d.f90wrap_crosssection__set__coord(self._handle, coord)
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 81
        
        """
        return _dassflow1d.f90wrap_crosssection__get__x(self._handle)
    
    @x.setter
    def x(self, x):
        _dassflow1d.f90wrap_crosssection__set__x(self._handle, x)
    
    @property
    def level(self):
        """
        Element level ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 83
        
        """
        return _dassflow1d.f90wrap_crosssection__get__level(self._handle)
    
    @level.setter
    def level(self, level):
        _dassflow1d.f90wrap_crosssection__set__level(self._handle, level)
    
    @property
    def nlevels(self):
        """
        Element nlevels ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 85
        
        """
        return _dassflow1d.f90wrap_crosssection__get__nlevels(self._handle)
    
    @nlevels.setter
    def nlevels(self, nlevels):
        _dassflow1d.f90wrap_crosssection__set__nlevels(self._handle, nlevels)
    
    @property
    def bathy(self):
        """
        Element bathy ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 91
        
        """
        return _dassflow1d.f90wrap_crosssection__get__bathy(self._handle)
    
    @bathy.setter
    def bathy(self, bathy):
        _dassflow1d.f90wrap_crosssection__set__bathy(self._handle, bathy)
    
    @property
    def delta(self):
        """
        Element delta ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 93
        
        """
        return _dassflow1d.f90wrap_crosssection__get__delta(self._handle)
    
    @delta.setter
    def delta(self, delta):
        _dassflow1d.f90wrap_crosssection__set__delta(self._handle, delta)
    
    @property
    def deltademi(self):
        """
        Element deltademi ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 95
        
        """
        return _dassflow1d.f90wrap_crosssection__get__deltademi(self._handle)
    
    @deltademi.setter
    def deltademi(self, deltademi):
        _dassflow1d.f90wrap_crosssection__set__deltademi(self._handle, deltademi)
    
    @property
    def slope(self):
        """
        Element slope ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 97
        
        """
        return _dassflow1d.f90wrap_crosssection__get__slope(self._handle)
    
    @slope.setter
    def slope(self, slope):
        _dassflow1d.f90wrap_crosssection__set__slope(self._handle, slope)
    
    @property
    def level_heights(self):
        """
        Element level_heights ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 99
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__level_heights(self._handle)
        if array_handle in self._arrays:
            level_heights = self._arrays[array_handle]
        else:
            level_heights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__level_heights)
            self._arrays[array_handle] = level_heights
        return level_heights
    
    @level_heights.setter
    def level_heights(self, level_heights):
        self.level_heights[...] = level_heights
    
    @property
    def level_widths(self):
        """
        Element level_widths ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 101
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__level_widths(self._handle)
        if array_handle in self._arrays:
            level_widths = self._arrays[array_handle]
        else:
            level_widths = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__level_widths)
            self._arrays[array_handle] = level_widths
        return level_widths
    
    @level_widths.setter
    def level_widths(self, level_widths):
        self.level_widths[...] = level_widths
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 103
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__y(self._handle)
        if array_handle in self._arrays:
            y = self._arrays[array_handle]
        else:
            y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__y)
            self._arrays[array_handle] = y
        return y
    
    @y.setter
    def y(self, y):
        self.y[...] = y
    
    @property
    def strickler_params(self):
        """
        Element strickler_params ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 105
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__strickler_params(self._handle)
        if array_handle in self._arrays:
            strickler_params = self._arrays[array_handle]
        else:
            strickler_params = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__strickler_params)
            self._arrays[array_handle] = strickler_params
        return strickler_params
    
    @strickler_params.setter
    def strickler_params(self, strickler_params):
        self.strickler_params[...] = strickler_params
    
    @property
    def poly(self):
        """
        Element poly ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 108
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__poly(self._handle)
        if array_handle in self._arrays:
            poly = self._arrays[array_handle]
        else:
            poly = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__poly)
            self._arrays[array_handle] = poly
        return poly
    
    @poly.setter
    def poly(self, poly):
        self.poly[...] = poly
    
    @property
    def area_cum(self):
        """
        Element area_cum ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 110
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__area_cum(self._handle)
        if array_handle in self._arrays:
            area_cum = self._arrays[array_handle]
        else:
            area_cum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__area_cum)
            self._arrays[array_handle] = area_cum
        return area_cum
    
    @area_cum.setter
    def area_cum(self, area_cum):
        self.area_cum[...] = area_cum
    
    @property
    def perim_cum(self):
        """
        Element perim_cum ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 112
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__perim_cum(self._handle)
        if array_handle in self._arrays:
            perim_cum = self._arrays[array_handle]
        else:
            perim_cum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__perim_cum)
            self._arrays[array_handle] = perim_cum
        return perim_cum
    
    @perim_cum.setter
    def perim_cum(self, perim_cum):
        self.perim_cum[...] = perim_cum
    
    @property
    def pa_cum(self):
        """
        Element pa_cum ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 114
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_crosssection__array__pa_cum(self._handle)
        if array_handle in self._arrays:
            pa_cum = self._arrays[array_handle]
        else:
            pa_cum = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_crosssection__array__pa_cum)
            self._arrays[array_handle] = pa_cum
        return pa_cum
    
    @pa_cum.setter
    def pa_cum(self, pa_cum):
        self.pa_cum[...] = pa_cum
    
    def __str__(self):
        ret = ['<crosssection>{\n']
        ret.append('    coord : ')
        ret.append(repr(self.coord))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    level : ')
        ret.append(repr(self.level))
        ret.append(',\n    nlevels : ')
        ret.append(repr(self.nlevels))
        ret.append(',\n    bathy : ')
        ret.append(repr(self.bathy))
        ret.append(',\n    delta : ')
        ret.append(repr(self.delta))
        ret.append(',\n    deltademi : ')
        ret.append(repr(self.deltademi))
        ret.append(',\n    slope : ')
        ret.append(repr(self.slope))
        ret.append(',\n    level_heights : ')
        ret.append(repr(self.level_heights))
        ret.append(',\n    level_widths : ')
        ret.append(repr(self.level_widths))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append(',\n    strickler_params : ')
        ret.append(repr(self.strickler_params))
        ret.append(',\n    poly : ')
        ret.append(repr(self.poly))
        ret.append(',\n    area_cum : ')
        ret.append(repr(self.area_cum))
        ret.append(',\n    perim_cum : ')
        ret.append(repr(self.perim_cum))
        ret.append(',\n    pa_cum : ')
        ret.append(repr(self.pa_cum))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Segment")
class Segment(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=segment)
    
    
    Defined at m_mesh.f90 lines 117-129
    
    """
    def __init__(self, nus_segs, handle=None):
        """
        self = Segment(nus_segs)
        
        
        Defined at m_mesh.f90 lines 559-565
        
        Parameters
        ----------
        nus_segs : int
        
        Returns
        -------
        seg : Segment
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_segment_initialise(nus_segs=nus_segs)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Segment
        
        
        Defined at m_mesh.f90 lines 567-570
        
        Parameters
        ----------
        seg : Segment
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_segment_finalise(seg=self._handle)
    
    @property
    def first_cs(self):
        """
        Element first_cs ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 119
        
        """
        return _dassflow1d.f90wrap_segment__get__first_cs(self._handle)
    
    @first_cs.setter
    def first_cs(self, first_cs):
        _dassflow1d.f90wrap_segment__set__first_cs(self._handle, first_cs)
    
    @property
    def last_cs(self):
        """
        Element last_cs ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 121
        
        """
        return _dassflow1d.f90wrap_segment__get__last_cs(self._handle)
    
    @last_cs.setter
    def last_cs(self, last_cs):
        _dassflow1d.f90wrap_segment__set__last_cs(self._handle, last_cs)
    
    @property
    def ds_seg(self):
        """
        Element ds_seg ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 123
        
        """
        return _dassflow1d.f90wrap_segment__get__ds_seg(self._handle)
    
    @ds_seg.setter
    def ds_seg(self, ds_seg):
        _dassflow1d.f90wrap_segment__set__ds_seg(self._handle, ds_seg)
    
    @property
    def us_seg(self):
        """
        Element us_seg ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 125
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_segment__array__us_seg(self._handle)
        if array_handle in self._arrays:
            us_seg = self._arrays[array_handle]
        else:
            us_seg = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_segment__array__us_seg)
            self._arrays[array_handle] = us_seg
        return us_seg
    
    @us_seg.setter
    def us_seg(self, us_seg):
        self.us_seg[...] = us_seg
    
    @property
    def ds_bc(self):
        """
        Element ds_bc ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 127
        
        """
        return _dassflow1d.f90wrap_segment__get__ds_bc(self._handle)
    
    @ds_bc.setter
    def ds_bc(self, ds_bc):
        _dassflow1d.f90wrap_segment__set__ds_bc(self._handle, ds_bc)
    
    @property
    def us_bc(self):
        """
        Element us_bc ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 129
        
        """
        return _dassflow1d.f90wrap_segment__get__us_bc(self._handle)
    
    @us_bc.setter
    def us_bc(self, us_bc):
        _dassflow1d.f90wrap_segment__set__us_bc(self._handle, us_bc)
    
    def __str__(self):
        ret = ['<segment>{\n']
        ret.append('    first_cs : ')
        ret.append(repr(self.first_cs))
        ret.append(',\n    last_cs : ')
        ret.append(repr(self.last_cs))
        ret.append(',\n    ds_seg : ')
        ret.append(repr(self.ds_seg))
        ret.append(',\n    us_seg : ')
        ret.append(repr(self.us_seg))
        ret.append(',\n    ds_bc : ')
        ret.append(repr(self.ds_bc))
        ret.append(',\n    us_bc : ')
        ret.append(repr(self.us_bc))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.SpatialField")
class SpatialField(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=spatialfield)
    
    
    Defined at m_mesh.f90 lines 132-135
    
    """
    def __init__(self, handle=None):
        """
        self = Spatialfield()
        
        
        Defined at m_mesh.f90 lines 132-135
        
        
        Returns
        -------
        this : Spatialfield
        	Object to be constructed
        
        
        Automatically generated constructor for spatialfield
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_spatialfield_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Spatialfield
        
        
        Defined at m_mesh.f90 lines 132-135
        
        Parameters
        ----------
        this : Spatialfield
        	Object to be destructed
        
        
        Automatically generated destructor for spatialfield
        """
        if self._alloc:
            _dassflow1d.f90wrap_spatialfield_finalise(this=self._handle)
    
    @property
    def interp(self):
        """
        Element interp ftype=character(len=8) pytype=str
        
        
        Defined at m_mesh.f90 line 133
        
        """
        return _dassflow1d.f90wrap_spatialfield__get__interp(self._handle)
    
    @interp.setter
    def interp(self, interp):
        _dassflow1d.f90wrap_spatialfield__set__interp(self._handle, interp)
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 134
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_spatialfield__array__x(self._handle)
        if array_handle in self._arrays:
            x = self._arrays[array_handle]
        else:
            x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_spatialfield__array__x)
            self._arrays[array_handle] = x
        return x
    
    @x.setter
    def x(self, x):
        self.x[...] = x
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 135
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_spatialfield__array__y(self._handle)
        if array_handle in self._arrays:
            y = self._arrays[array_handle]
        else:
            y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_spatialfield__array__y)
            self._arrays[array_handle] = y
        return y
    
    @y.setter
    def y(self, y):
        self.y[...] = y
    
    def __str__(self):
        ret = ['<spatialfield>{\n']
        ret.append('    interp : ')
        ret.append(repr(self.interp))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Mesh")
class Mesh(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=mesh)
    
    
    Defined at m_mesh.f90 lines 138-154
    
    """
    def dealloc_mesh(self):
        """
        dealloc_mesh(self)
        
        
        Defined at m_mesh.f90 lines 581-593
        
        Parameters
        ----------
        msh : Mesh
        
        """
        _dassflow1d.f90wrap_dealloc_mesh(msh=self._handle)
    
    def __init__(self, ncs, nseg=None, handle=None):
        """
        self = Mesh(ncs[, nseg])
        
        
        Defined at m_mesh.f90 lines 606-616
        
        Parameters
        ----------
        ncs : int
        nseg : int
        
        Returns
        -------
        msh : Mesh
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_mesh_initialise(ncs=ncs, nseg=nseg)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def setup_segment(self, iseg, first_cs, ncs, us_seg, ds_seg):
        """
        setup_segment(self, iseg, first_cs, ncs, us_seg, ds_seg)
        
        
        Defined at m_mesh.f90 lines 618-650
        
        Parameters
        ----------
        msh : Mesh
        iseg : int
        first_cs : int
        ncs : int
        us_seg : int array
        ds_seg : int
        
        """
        _dassflow1d.f90wrap_setup_segment(msh=self._handle, iseg=iseg, \
            first_cs=first_cs, ncs=ncs, us_seg=us_seg, ds_seg=ds_seg)
    
    def setup_crosssection(self, ics, nlevels):
        """
        setup_crosssection(self, ics, nlevels)
        
        
        Defined at m_mesh.f90 lines 652-657
        
        Parameters
        ----------
        msh : Mesh
        ics : int
        nlevels : int
        
        """
        _dassflow1d.f90wrap_setup_crosssection(msh=self._handle, ics=ics, \
            nlevels=nlevels)
    
    def update_geometries(self):
        """
        update_geometries(self)
        
        
        Defined at m_mesh.f90 lines 659-666
        
        Parameters
        ----------
        msh : Mesh
        
        """
        _dassflow1d.f90wrap_update_geometries(msh=self._handle)
    
    def set_strickler_type(self, strickler_type):
        """
        set_strickler_type(self, strickler_type)
        
        
        Defined at m_mesh.f90 lines 668-679
        
        Parameters
        ----------
        msh : Mesh
        strickler_type : str
        
        """
        _dassflow1d.f90wrap_set_strickler_type(msh=self._handle, \
            strickler_type=strickler_type)
    
    def set_uniform_strickler_parameters(self, strickler_params):
        """
        set_uniform_strickler_parameters(self, strickler_params)
        
        
        Defined at m_mesh.f90 lines 681-689
        
        Parameters
        ----------
        msh : Mesh
        strickler_params : float array
        
        """
        _dassflow1d.f90wrap_set_uniform_strickler_parameters(msh=self._handle, \
            strickler_params=strickler_params)
    
    def set_bathy_field_linear(self, x, bathy):
        """
        set_bathy_field_linear(self, x, bathy)
        
        
        Defined at m_mesh.f90 lines 691-705
        
        Parameters
        ----------
        msh : Mesh
        x : float array
        bathy : float array
        
        """
        _dassflow1d.f90wrap_set_bathy_field_linear(msh=self._handle, x=x, bathy=bathy)
    
    def set_strickler_fields_segment(self, strickler_params):
        """
        set_strickler_fields_segment(self, strickler_params)
        
        
        Defined at m_mesh.f90 lines 707-755
        
        Parameters
        ----------
        msh : Mesh
        strickler_params : float array
        
        """
        _dassflow1d.f90wrap_set_strickler_fields_segment(msh=self._handle, \
            strickler_params=strickler_params)
    
    def set_strickler_fields_linear(self, x, strickler_params):
        """
        set_strickler_fields_linear(self, x, strickler_params)
        
        
        Defined at m_mesh.f90 lines 757-804
        
        Parameters
        ----------
        msh : Mesh
        x : float array
        strickler_params : float array
        
        """
        _dassflow1d.f90wrap_set_strickler_fields_linear(msh=self._handle, x=x, \
            strickler_params=strickler_params)
    
    def get_segment_cs_index(self, iseg, ics_seg):
        """
        ics = get_segment_cs_index(self, iseg, ics_seg)
        
        
        Defined at m_mesh.f90 lines 820-826
        
        Parameters
        ----------
        msh : Mesh
        iseg : int
        ics_seg : int
        
        Returns
        -------
        ics : int
        
        """
        ics = _dassflow1d.f90wrap_get_segment_cs_index(msh=self._handle, iseg=iseg, \
            ics_seg=ics_seg)
        return ics
    
    def __del__(self):
        """
        Destructor for class Mesh
        
        
        Defined at m_mesh.f90 lines 138-154
        
        Parameters
        ----------
        this : Mesh
        	Object to be destructed
        
        
        Automatically generated destructor for mesh
        """
        if self._alloc:
            _dassflow1d.f90wrap_mesh_finalise(this=self._handle)
    
    @property
    def ncs(self):
        """
        Element ncs ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 140
        
        """
        return _dassflow1d.f90wrap_mesh__get__ncs(self._handle)
    
    @ncs.setter
    def ncs(self, ncs):
        _dassflow1d.f90wrap_mesh__set__ncs(self._handle, ncs)
    
    @property
    def nseg(self):
        """
        Element nseg ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 142
        
        """
        return _dassflow1d.f90wrap_mesh__get__nseg(self._handle)
    
    @nseg.setter
    def nseg(self, nseg):
        _dassflow1d.f90wrap_mesh__set__nseg(self._handle, nseg)
    
    def init_array_cs(self):
        self.cs = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_mesh__array_getitem__cs,
                                        _dassflow1d.f90wrap_mesh__array_setitem__cs,
                                        _dassflow1d.f90wrap_mesh__array_len__cs,
                                        """
        Element cs ftype=type(crosssection) pytype=Crosssection
        
        
        Defined at m_mesh.f90 line 144
        
        """, Crosssection)
        return self.cs
    
    def init_array_seg(self):
        self.seg = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_mesh__array_getitem__seg,
                                        _dassflow1d.f90wrap_mesh__array_setitem__seg,
                                        _dassflow1d.f90wrap_mesh__array_len__seg,
                                        """
        Element seg ftype=type(segment) pytype=Segment
        
        
        Defined at m_mesh.f90 line 146
        
        """, Segment)
        return self.seg
    
    @property
    def bathy_field(self):
        """
        Element bathy_field ftype=type(spatialfield) pytype=Spatialfield
        
        
        Defined at m_mesh.f90 line 147
        
        """
        bathy_field_handle = _dassflow1d.f90wrap_mesh__get__bathy_field(self._handle)
        if tuple(bathy_field_handle) in self._objs:
            bathy_field = self._objs[tuple(bathy_field_handle)]
        else:
            bathy_field = SpatialField.from_handle(bathy_field_handle)
            self._objs[tuple(bathy_field_handle)] = bathy_field
        return bathy_field
    
    @bathy_field.setter
    def bathy_field(self, bathy_field):
        bathy_field = bathy_field._handle
        _dassflow1d.f90wrap_mesh__set__bathy_field(self._handle, bathy_field)
    
    @property
    def strickler_type(self):
        """
        Element strickler_type ftype=character(len=16) pytype=str
        
        
        Defined at m_mesh.f90 line 148
        
        """
        return _dassflow1d.f90wrap_mesh__get__strickler_type(self._handle)
    
    @strickler_type.setter
    def strickler_type(self, strickler_type):
        _dassflow1d.f90wrap_mesh__set__strickler_type(self._handle, strickler_type)
    
    @property
    def strickler_type_code(self):
        """
        Element strickler_type_code ftype=integer  pytype=int
        
        
        Defined at m_mesh.f90 line 149
        
        """
        return _dassflow1d.f90wrap_mesh__get__strickler_type_code(self._handle)
    
    @strickler_type_code.setter
    def strickler_type_code(self, strickler_type_code):
        _dassflow1d.f90wrap_mesh__set__strickler_type_code(self._handle, \
            strickler_type_code)
    
    def init_array_strickler_fields(self):
        self.strickler_fields = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_mesh__array_getitem__strickler_fields,
                                        _dassflow1d.f90wrap_mesh__array_setitem__strickler_fields,
                                        _dassflow1d.f90wrap_mesh__array_len__strickler_fields,
                                        """
        Element strickler_fields ftype=type(spatialfield) pytype=Spatialfield
        
        
        Defined at m_mesh.f90 line 150
        
        """, SpatialField)
        return self.strickler_fields
    
    def __str__(self):
        ret = ['<mesh>{\n']
        ret.append('    ncs : ')
        ret.append(repr(self.ncs))
        ret.append(',\n    nseg : ')
        ret.append(repr(self.nseg))
        ret.append(',\n    bathy_field : ')
        ret.append(repr(self.bathy_field))
        ret.append(',\n    strickler_type : ')
        ret.append(repr(self.strickler_type))
        ret.append(',\n    strickler_type_code : ')
        ret.append(repr(self.strickler_type_code))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_cs, init_array_seg, \
        init_array_strickler_fields]
    

@f90wrap.runtime.register_class("dassflow1d.point_in_mesh")
class point_in_mesh(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=point_in_mesh)
    
    
    Defined at m_mesh.f90 lines 159-161
    
    """
    def __init__(self, handle=None):
        """
        self = Point_In_Mesh()
        
        
        Defined at m_mesh.f90 lines 159-161
        
        
        Returns
        -------
        this : Point_In_Mesh
        	Object to be constructed
        
        
        Automatically generated constructor for point_in_mesh
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_point_in_mesh_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Point_In_Mesh
        
        
        Defined at m_mesh.f90 lines 159-161
        
        Parameters
        ----------
        this : Point_In_Mesh
        	Object to be destructed
        
        
        Automatically generated destructor for point_in_mesh
        """
        if self._alloc:
            _dassflow1d.f90wrap_point_in_mesh_finalise(this=self._handle)
    
    @property
    def indexi(self):
        """
        Element indexi ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 161
        
        """
        return _dassflow1d.f90wrap_point_in_mesh__get__indexi(self._handle)
    
    @indexi.setter
    def indexi(self, indexi):
        _dassflow1d.f90wrap_point_in_mesh__set__indexi(self._handle, indexi)
    
    def __str__(self):
        ret = ['<point_in_mesh>{\n']
        ret.append('    indexi : ')
        ret.append(repr(self.indexi))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def get_k_constant():
    """
    Element k_constant ftype=integer pytype=int
    
    
    Defined at m_mesh.f90 line 72
    
    """
    return _dassflow1d.f90wrap_m_mesh__get__k_constant()

K_constant = get_k_constant()

def get_k_powerlaw_h():
    """
    Element k_powerlaw_h ftype=integer pytype=int
    
    
    Defined at m_mesh.f90 line 73
    
    """
    return _dassflow1d.f90wrap_m_mesh__get__k_powerlaw_h()

K_powerlaw_h = get_k_powerlaw_h()

def get_strickler_type_constant():
    """
    Element strickler_type_constant ftype=integer pytype=int
    
    
    Defined at m_mesh.f90 line 74
    
    """
    return _dassflow1d.f90wrap_m_mesh__get__strickler_type_constant()

strickler_type_constant = get_strickler_type_constant()

def get_strickler_type_powerlaw_h():
    """
    Element strickler_type_powerlaw_h ftype=integer pytype=int
    
    
    Defined at m_mesh.f90 line 75
    
    """
    return _dassflow1d.f90wrap_m_mesh__get__strickler_type_powerlaw_h()

strickler_type_powerlaw_h = get_strickler_type_powerlaw_h()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_mesh".')

for func in _dt_array_initialisers:
    func()
