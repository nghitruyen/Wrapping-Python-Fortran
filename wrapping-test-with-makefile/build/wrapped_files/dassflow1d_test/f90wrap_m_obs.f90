! Module m_obs defined in file m_obs.f90

subroutine f90wrap_obsstation__array__ics(this, nd, dtype, dshape, dloc)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ics)) then
        dshape(1:1) = shape(this_ptr%p%ics)
        dloc = loc(this_ptr%p%ics)
    else
        dloc = 0
    end if
end subroutine f90wrap_obsstation__array__ics

subroutine f90wrap_obsstation__get__iobs(this, f90wrap_iobs)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in)   :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_iobs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_iobs = this_ptr%p%iobs
end subroutine f90wrap_obsstation__get__iobs

subroutine f90wrap_obsstation__set__iobs(this, f90wrap_iobs)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in)   :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_iobs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%iobs = f90wrap_iobs
end subroutine f90wrap_obsstation__set__iobs

subroutine f90wrap_obsstation__array__t(this, nd, dtype, dshape, dloc)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_obsstation__array__t

subroutine f90wrap_obsstation__get__offset(this, f90wrap_offset)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in)   :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_offset
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_offset = this_ptr%p%offset
end subroutine f90wrap_obsstation__get__offset

subroutine f90wrap_obsstation__set__offset(this, f90wrap_offset)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in)   :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_offset
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%offset = f90wrap_offset
end subroutine f90wrap_obsstation__set__offset

subroutine f90wrap_obsstation__array__obs(this, nd, dtype, dshape, dloc)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%obs)) then
        dshape(1:2) = shape(this_ptr%p%obs)
        dloc = loc(this_ptr%p%obs)
    else
        dloc = 0
    end if
end subroutine f90wrap_obsstation__array__obs

subroutine f90wrap_obsstation__array__w(this, nd, dtype, dshape, dloc)
    use m_obs, only: obsstation
    implicit none
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: this(2)
    type(obsstation_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%w)) then
        dshape(1:2) = shape(this_ptr%p%w)
        dloc = loc(this_ptr%p%w)
    else
        dloc = 0
    end if
end subroutine f90wrap_obsstation__array__w

subroutine f90wrap_obsstation_initialise(station, ics, t, values, n0, n1, n2, n3)
    use m_obs, only: obsstation_initialise, obsstation
    implicit none
    
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    type(obsstation_ptr_type) :: station_ptr
    integer, intent(out), dimension(2) :: station
    integer(4), intent(in), dimension(n0) :: ics
    real(8), intent(in), dimension(n1) :: t
    real(8), intent(in), dimension(n2,n3) :: values
    integer :: n0
    !f2py intent(hide), depend(ics) :: n0 = shape(ics,0)
    integer :: n1
    !f2py intent(hide), depend(t) :: n1 = shape(t,0)
    integer :: n2
    !f2py intent(hide), depend(values) :: n2 = shape(values,0)
    integer :: n3
    !f2py intent(hide), depend(values) :: n3 = shape(values,1)
    allocate(station_ptr%p)
    call obsstation_initialise(station=station_ptr%p, ics=ics, t=t, values=values)
    station = transfer(station_ptr, station)
end subroutine f90wrap_obsstation_initialise

subroutine f90wrap_obsstation_finalise(station)
    use m_obs, only: obsstation, obsstation_finalise
    implicit none
    
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    type(obsstation_ptr_type) :: station_ptr
    integer, intent(in), dimension(2) :: station
    station_ptr = transfer(station, station_ptr)
    call obsstation_finalise(station=station_ptr%p)
    deallocate(station_ptr%p)
end subroutine f90wrap_obsstation_finalise

subroutine f90wrap_setup_station(station, ics, t, values, n0, n1, n2, n3)
    use m_obs, only: obsstation, setup_station
    implicit none
    
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    type(obsstation_ptr_type) :: station_ptr
    integer, intent(in), dimension(2) :: station
    integer(4), intent(in), dimension(n0) :: ics
    real(8), intent(in), dimension(n1) :: t
    real(8), intent(in), optional, dimension(n2,n3) :: values
    integer :: n0
    !f2py intent(hide), depend(ics) :: n0 = shape(ics,0)
    integer :: n1
    !f2py intent(hide), depend(t) :: n1 = shape(t,0)
    integer :: n2
    !f2py intent(hide), depend(values) :: n2 = shape(values,0)
    integer :: n3
    !f2py intent(hide), depend(values) :: n3 = shape(values,1)
    station_ptr = transfer(station, station_ptr)
    call setup_station(station=station_ptr%p, ics=ics, t=t, values=values)
end subroutine f90wrap_setup_station

subroutine f90wrap_set_station_weights(station, weights, n0, n1)
    use m_obs, only: obsstation, set_station_weights
    implicit none
    
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    type(obsstation_ptr_type) :: station_ptr
    integer, intent(in), dimension(2) :: station
    real(8), intent(in), dimension(n0,n1) :: weights
    integer :: n0
    !f2py intent(hide), depend(weights) :: n0 = shape(weights,0)
    integer :: n1
    !f2py intent(hide), depend(weights) :: n1 = shape(weights,1)
    station_ptr = transfer(station, station_ptr)
    call set_station_weights(station=station_ptr%p, weights=weights)
end subroutine f90wrap_set_station_weights

subroutine f90wrap_observations__array_getitem__stations(f90wrap_this, f90wrap_i, stationsitem)
    
    use m_obs, only: obsstation, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(observations_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: stationsitem(2)
    type(obsstation_ptr_type) :: stations_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%stations)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%stations)) then
            call f90wrap_abort("array index out of range")
        else
            stations_ptr%p => this_ptr%p%stations(f90wrap_i)
            stationsitem = transfer(stations_ptr,stationsitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_observations__array_getitem__stations

subroutine f90wrap_observations__array_setitem__stations(f90wrap_this, f90wrap_i, stationsitem)
    
    use m_obs, only: obsstation, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(observations_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: stationsitem(2)
    type(obsstation_ptr_type) :: stations_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%stations)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%stations)) then
            call f90wrap_abort("array index out of range")
        else
            stations_ptr = transfer(stationsitem,stations_ptr)
            this_ptr%p%stations(f90wrap_i) = stations_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_observations__array_setitem__stations

subroutine f90wrap_observations__array_len__stations(f90wrap_this, f90wrap_n)
    
    use m_obs, only: obsstation, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type obsstation_ptr_type
        type(obsstation), pointer :: p => NULL()
    end type obsstation_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(observations_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%stations)) then
        f90wrap_n = size(this_ptr%p%stations)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_observations__array_len__stations

subroutine f90wrap_observations__array__obs(this, nd, dtype, dshape, dloc)
    use m_obs, only: observations
    implicit none
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    integer, intent(in) :: this(2)
    type(observations_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%obs)) then
        dshape(1:2) = shape(this_ptr%p%obs)
        dloc = loc(this_ptr%p%obs)
    else
        dloc = 0
    end if
end subroutine f90wrap_observations__array__obs

subroutine f90wrap_observations__array__w(this, nd, dtype, dshape, dloc)
    use m_obs, only: observations
    implicit none
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    integer, intent(in) :: this(2)
    type(observations_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%w)) then
        dshape(1:2) = shape(this_ptr%p%w)
        dloc = loc(this_ptr%p%w)
    else
        dloc = 0
    end if
end subroutine f90wrap_observations__array__w

subroutine f90wrap_observations__array__est(this, nd, dtype, dshape, dloc)
    use m_obs, only: observations
    implicit none
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    integer, intent(in) :: this(2)
    type(observations_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%est)) then
        dshape(1:2) = shape(this_ptr%p%est)
        dloc = loc(this_ptr%p%est)
    else
        dloc = 0
    end if
end subroutine f90wrap_observations__array__est

subroutine f90wrap_observations_initialise(obs, nstations)
    use m_obs, only: observations, observations_initialise
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type(observations_ptr_type) :: obs_ptr
    integer, intent(out), dimension(2) :: obs
    integer(4), intent(in) :: nstations
    allocate(obs_ptr%p)
    call observations_initialise(obs=obs_ptr%p, nstations=nstations)
    obs = transfer(obs_ptr, obs)
end subroutine f90wrap_observations_initialise

subroutine f90wrap_observations_finalise(obs)
    use m_obs, only: observations_finalise, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type(observations_ptr_type) :: obs_ptr
    integer, intent(in), dimension(2) :: obs
    obs_ptr = transfer(obs, obs_ptr)
    call observations_finalise(obs=obs_ptr%p)
    deallocate(obs_ptr%p)
end subroutine f90wrap_observations_finalise

subroutine f90wrap_setup_observations_data(obs, ndata)
    use m_obs, only: setup_observations_data, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type(observations_ptr_type) :: obs_ptr
    integer, intent(in), dimension(2) :: obs
    integer(4), intent(in) :: ndata
    obs_ptr = transfer(obs, obs_ptr)
    call setup_observations_data(obs=obs_ptr%p, ndata=ndata)
end subroutine f90wrap_setup_observations_data

subroutine f90wrap_data_from_stations(obs, discharge_estimation)
    use m_obs, only: data_from_stations, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type(observations_ptr_type) :: obs_ptr
    integer, intent(in), dimension(2) :: obs
    logical, intent(in), optional :: discharge_estimation
    obs_ptr = transfer(obs, obs_ptr)
    call data_from_stations(obs=obs_ptr%p, discharge_estimation=discharge_estimation)
end subroutine f90wrap_data_from_stations

subroutine f90wrap_reset_observations_counters(obs)
    use m_obs, only: reset_observations_counters, observations
    implicit none
    
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    type(observations_ptr_type) :: obs_ptr
    integer, intent(in), dimension(2) :: obs
    obs_ptr = transfer(obs, obs_ptr)
    call reset_observations_counters(obs=obs_ptr%p)
end subroutine f90wrap_reset_observations_counters

subroutine f90wrap_all_observed(t, msh, obs, n0)
    use m_obs, only: all_observed, observations
    use m_mesh, only: mesh
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type observations_ptr_type
        type(observations), pointer :: p => NULL()
    end type observations_ptr_type
    real(8), intent(in), dimension(n0) :: t
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    type(observations_ptr_type) :: obs_ptr
    integer, intent(out), dimension(2) :: obs
    integer :: n0
    !f2py intent(hide), depend(t) :: n0 = shape(t,0)
    msh_ptr = transfer(msh, msh_ptr)
    allocate(obs_ptr%p)
    call all_observed(t=t, msh=msh_ptr%p, obs=obs_ptr%p)
    obs = transfer(obs_ptr, obs)
end subroutine f90wrap_all_observed

subroutine f90wrap_cost_l2_heights(hest, ret_cost, hobs, n0, n1)
    use m_obs, only: cost_l2_heights
    implicit none
    
    real(8), intent(in), dimension(n0) :: hest
    real(8), intent(out) :: ret_cost
    real(8), intent(in), dimension(n1) :: hobs
    integer :: n0
    !f2py intent(hide), depend(hest) :: n0 = shape(hest,0)
    integer :: n1
    !f2py intent(hide), depend(hobs) :: n1 = shape(hobs,0)
    ret_cost = cost_l2_heights(Hest=hest, Hobs=hobs)
end subroutine f90wrap_cost_l2_heights

subroutine f90wrap_cost_r2_heights(hest, hobs, ret_cost, wobs, n0, n1, n2)
    use m_obs, only: cost_r2_heights
    implicit none
    
    real(8), intent(in), dimension(n0) :: hest
    real(8), intent(in), dimension(n1) :: hobs
    real(8), intent(out) :: ret_cost
    real(8), intent(in), dimension(n2) :: wobs
    integer :: n0
    !f2py intent(hide), depend(hest) :: n0 = shape(hest,0)
    integer :: n1
    !f2py intent(hide), depend(hobs) :: n1 = shape(hobs,0)
    integer :: n2
    !f2py intent(hide), depend(wobs) :: n2 = shape(wobs,0)
    ret_cost = cost_r2_heights(Hest=hest, Hobs=hobs, wobs=wobs)
end subroutine f90wrap_cost_r2_heights

! End of module m_obs defined in file m_obs.f90

