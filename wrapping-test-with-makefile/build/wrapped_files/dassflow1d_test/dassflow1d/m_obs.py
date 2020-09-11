"""
Module m_obs


Defined at m_obs.f90 lines 68-408

"""
from __future__ import print_function, absolute_import, division
import _dassflow1d
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("dassflow1d.ObsStation")
class ObsStation(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=obsstation)
    
    
    Defined at m_obs.f90 lines 72-84
    
    """
    def __init__(self, ics, t, values, handle=None):
        """
        self = Obsstation(ics, t, values)
        
        
        Defined at m_obs.f90 lines 97-116
        
        Parameters
        ----------
        ics : int array
        t : float array
        values : float array
        
        Returns
        -------
        station : Obsstation
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_obsstation_initialise(ics=ics, t=t, values=values)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Obsstation
        
        
        Defined at m_obs.f90 lines 118-124
        
        Parameters
        ----------
        station : Obsstation
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_obsstation_finalise(station=self._handle)
    
    def setup_station(self, ics, t, values=None):
        """
        setup_station(self, ics, t[, values])
        
        
        Defined at m_obs.f90 lines 126-156
        
        Parameters
        ----------
        station : Obsstation
        ics : int array
        t : float array
        values : float array
        
        """
        _dassflow1d.f90wrap_setup_station(station=self._handle, ics=ics, t=t, \
            values=values)
    
    def set_station_weights(self, weights):
        """
        set_station_weights(self, weights)
        
        
        Defined at m_obs.f90 lines 158-171
        
        Parameters
        ----------
        station : Obsstation
        weights : float array
        
        """
        _dassflow1d.f90wrap_set_station_weights(station=self._handle, weights=weights)
    
    @property
    def ics(self):
        """
        Element ics ftype=integer(ip) pytype=int
        
        
        Defined at m_obs.f90 line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_obsstation__array__ics(self._handle)
        if array_handle in self._arrays:
            ics = self._arrays[array_handle]
        else:
            ics = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_obsstation__array__ics)
            self._arrays[array_handle] = ics
        return ics
    
    @ics.setter
    def ics(self, ics):
        self.ics[...] = ics
    
    @property
    def iobs(self):
        """
        Element iobs ftype=integer(ip) pytype=int
        
        
        Defined at m_obs.f90 line 76
        
        """
        return _dassflow1d.f90wrap_obsstation__get__iobs(self._handle)
    
    @iobs.setter
    def iobs(self, iobs):
        _dassflow1d.f90wrap_obsstation__set__iobs(self._handle, iobs)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 78
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_obsstation__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_obsstation__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def offset(self):
        """
        Element offset ftype=integer(ip) pytype=int
        
        
        Defined at m_obs.f90 line 80
        
        """
        return _dassflow1d.f90wrap_obsstation__get__offset(self._handle)
    
    @offset.setter
    def offset(self, offset):
        _dassflow1d.f90wrap_obsstation__set__offset(self._handle, offset)
    
    @property
    def obs(self):
        """
        Element obs ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 82
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_obsstation__array__obs(self._handle)
        if array_handle in self._arrays:
            obs = self._arrays[array_handle]
        else:
            obs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_obsstation__array__obs)
            self._arrays[array_handle] = obs
        return obs
    
    @obs.setter
    def obs(self, obs):
        self.obs[...] = obs
    
    @property
    def w(self):
        """
        Element w ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 84
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_obsstation__array__w(self._handle)
        if array_handle in self._arrays:
            w = self._arrays[array_handle]
        else:
            w = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_obsstation__array__w)
            self._arrays[array_handle] = w
        return w
    
    @w.setter
    def w(self, w):
        self.w[...] = w
    
    def __str__(self):
        ret = ['<obsstation>{\n']
        ret.append('    ics : ')
        ret.append(repr(self.ics))
        ret.append(',\n    iobs : ')
        ret.append(repr(self.iobs))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    offset : ')
        ret.append(repr(self.offset))
        ret.append(',\n    obs : ')
        ret.append(repr(self.obs))
        ret.append(',\n    w : ')
        ret.append(repr(self.w))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("dassflow1d.Observations")
class Observations(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=observations)
    
    
    Defined at m_obs.f90 lines 86-93
    
    """
    def __init__(self, nstations, handle=None):
        """
        self = Observations(nstations)
        
        
        Defined at m_obs.f90 lines 173-177
        
        Parameters
        ----------
        nstations : int
        
        Returns
        -------
        obs : Observations
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _dassflow1d.f90wrap_observations_initialise(nstations=nstations)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Observations
        
        
        Defined at m_obs.f90 lines 179-185
        
        Parameters
        ----------
        obs : Observations
        
        """
        if self._alloc:
            _dassflow1d.f90wrap_observations_finalise(obs=self._handle)
    
    def setup_observations_data(self, ndata):
        """
        setup_observations_data(self, ndata)
        
        
        Defined at m_obs.f90 lines 187-196
        
        Parameters
        ----------
        obs : Observations
        ndata : int
        
        """
        _dassflow1d.f90wrap_setup_observations_data(obs=self._handle, ndata=ndata)
    
    def data_from_stations(self, discharge_estimation=None):
        """
        data_from_stations(self[, discharge_estimation])
        
        
        Defined at m_obs.f90 lines 246-301
        
        Parameters
        ----------
        obs : Observations
        discharge_estimation : bool
        
        """
        _dassflow1d.f90wrap_data_from_stations(obs=self._handle, \
            discharge_estimation=discharge_estimation)
    
    def reset_observations_counters(self):
        """
        reset_observations_counters(self)
        
        
        Defined at m_obs.f90 lines 370-377
        
        Parameters
        ----------
        obs : Observations
        
        """
        _dassflow1d.f90wrap_reset_observations_counters(obs=self._handle)
    
    def init_array_stations(self):
        self.stations = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _dassflow1d.f90wrap_observations__array_getitem__stations,
                                        _dassflow1d.f90wrap_observations__array_setitem__stations,
                                        _dassflow1d.f90wrap_observations__array_len__stations,
                                        """
        Element stations ftype=type(obsstation) pytype=Obsstation
        
        
        Defined at m_obs.f90 line 87
        
        """, ObsStation)
        return self.stations
    
    @property
    def obs(self):
        """
        Element obs ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 89
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_observations__array__obs(self._handle)
        if array_handle in self._arrays:
            obs = self._arrays[array_handle]
        else:
            obs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_observations__array__obs)
            self._arrays[array_handle] = obs
        return obs
    
    @obs.setter
    def obs(self, obs):
        self.obs[...] = obs
    
    @property
    def w(self):
        """
        Element w ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 91
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_observations__array__w(self._handle)
        if array_handle in self._arrays:
            w = self._arrays[array_handle]
        else:
            w = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_observations__array__w)
            self._arrays[array_handle] = w
        return w
    
    @w.setter
    def w(self, w):
        self.w[...] = w
    
    @property
    def est(self):
        """
        Element est ftype=real(rp) pytype=float
        
        
        Defined at m_obs.f90 line 93
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _dassflow1d.f90wrap_observations__array__est(self._handle)
        if array_handle in self._arrays:
            est = self._arrays[array_handle]
        else:
            est = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _dassflow1d.f90wrap_observations__array__est)
            self._arrays[array_handle] = est
        return est
    
    @est.setter
    def est(self, est):
        self.est[...] = est
    
    def __str__(self):
        ret = ['<observations>{\n']
        ret.append('    obs : ')
        ret.append(repr(self.obs))
        ret.append(',\n    w : ')
        ret.append(repr(self.w))
        ret.append(',\n    est : ')
        ret.append(repr(self.est))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_stations]
    

def all_observed(t, msh):
    """
    obs = all_observed(t, msh)
    
    
    Defined at m_obs.f90 lines 198-244
    
    Parameters
    ----------
    t : float array
    msh : Mesh
    
    Returns
    -------
    obs : Observations
    
    """
    obs = _dassflow1d.f90wrap_all_observed(t=t, msh=msh._handle)
    obs = f90wrap.runtime.lookup_class("dassflow1d.Observations").from_handle(obs)
    return obs

def cost_l2_heights(hest, hobs):
    """
    cost = cost_l2_heights(hest, hobs)
    
    
    Defined at m_obs.f90 lines 379-392
    
    Parameters
    ----------
    hest : float array
    hobs : float array
    
    Returns
    -------
    cost : float
    
    """
    cost = _dassflow1d.f90wrap_cost_l2_heights(hest=hest, hobs=hobs)
    return cost

def cost_r2_heights(hest, hobs, wobs):
    """
    cost = cost_r2_heights(hest, hobs, wobs)
    
    
    Defined at m_obs.f90 lines 394-408
    
    Parameters
    ----------
    hest : float array
    hobs : float array
    wobs : float array
    
    Returns
    -------
    cost : float
    
    """
    cost = _dassflow1d.f90wrap_cost_r2_heights(hest=hest, hobs=hobs, wobs=wobs)
    return cost


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_obs".')

for func in _dt_array_initialisers:
    func()
