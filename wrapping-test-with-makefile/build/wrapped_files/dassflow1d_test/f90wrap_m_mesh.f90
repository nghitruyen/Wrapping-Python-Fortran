! Module m_mesh defined in file m_mesh.f90

subroutine f90wrap_crosssection__get__coord(this, f90wrap_coord)
    use m_mesh, only: crosssection
    use m_linear_algebra, only: vec2d
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr%p => this_ptr%p%coord
    f90wrap_coord = transfer(coord_ptr,f90wrap_coord)
end subroutine f90wrap_crosssection__get__coord

subroutine f90wrap_crosssection__set__coord(this, f90wrap_coord)
    use m_mesh, only: crosssection
    use m_linear_algebra, only: vec2d
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr = transfer(f90wrap_coord,coord_ptr)
    this_ptr%p%coord = coord_ptr%p
end subroutine f90wrap_crosssection__set__coord

subroutine f90wrap_crosssection__get__x(this, f90wrap_x)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_crosssection__get__x

subroutine f90wrap_crosssection__set__x(this, f90wrap_x)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x = f90wrap_x
end subroutine f90wrap_crosssection__set__x

subroutine f90wrap_crosssection__get__level(this, f90wrap_level)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_level
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_level = this_ptr%p%level
end subroutine f90wrap_crosssection__get__level

subroutine f90wrap_crosssection__set__level(this, f90wrap_level)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_level
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%level = f90wrap_level
end subroutine f90wrap_crosssection__set__level

subroutine f90wrap_crosssection__get__nlevels(this, f90wrap_nlevels)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nlevels
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nlevels = this_ptr%p%nlevels
end subroutine f90wrap_crosssection__get__nlevels

subroutine f90wrap_crosssection__set__nlevels(this, f90wrap_nlevels)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nlevels
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nlevels = f90wrap_nlevels
end subroutine f90wrap_crosssection__set__nlevels

subroutine f90wrap_crosssection__get__bathy(this, f90wrap_bathy)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_bathy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bathy = this_ptr%p%bathy
end subroutine f90wrap_crosssection__get__bathy

subroutine f90wrap_crosssection__set__bathy(this, f90wrap_bathy)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_bathy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bathy = f90wrap_bathy
end subroutine f90wrap_crosssection__set__bathy

subroutine f90wrap_crosssection__get__delta(this, f90wrap_delta)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_delta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_delta = this_ptr%p%delta
end subroutine f90wrap_crosssection__get__delta

subroutine f90wrap_crosssection__set__delta(this, f90wrap_delta)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_delta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%delta = f90wrap_delta
end subroutine f90wrap_crosssection__set__delta

subroutine f90wrap_crosssection__get__deltademi(this, f90wrap_deltademi)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_deltademi
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_deltademi = this_ptr%p%deltademi
end subroutine f90wrap_crosssection__get__deltademi

subroutine f90wrap_crosssection__set__deltademi(this, f90wrap_deltademi)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_deltademi
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%deltademi = f90wrap_deltademi
end subroutine f90wrap_crosssection__set__deltademi

subroutine f90wrap_crosssection__get__slope(this, f90wrap_slope)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_slope
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_slope = this_ptr%p%slope
end subroutine f90wrap_crosssection__get__slope

subroutine f90wrap_crosssection__set__slope(this, f90wrap_slope)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in)   :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_slope
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%slope = f90wrap_slope
end subroutine f90wrap_crosssection__set__slope

subroutine f90wrap_crosssection__array__level_heights(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%level_heights)) then
        dshape(1:1) = shape(this_ptr%p%level_heights)
        dloc = loc(this_ptr%p%level_heights)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__level_heights

subroutine f90wrap_crosssection__array__level_widths(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%level_widths)) then
        dshape(1:1) = shape(this_ptr%p%level_widths)
        dloc = loc(this_ptr%p%level_widths)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__level_widths

subroutine f90wrap_crosssection__array__y(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y)) then
        dshape(1:1) = shape(this_ptr%p%y)
        dloc = loc(this_ptr%p%y)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__y

subroutine f90wrap_crosssection__array__strickler_params(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%strickler_params)) then
        dshape(1:1) = shape(this_ptr%p%strickler_params)
        dloc = loc(this_ptr%p%strickler_params)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__strickler_params

subroutine f90wrap_crosssection__array__poly(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%poly)) then
        dshape(1:2) = shape(this_ptr%p%poly)
        dloc = loc(this_ptr%p%poly)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__poly

subroutine f90wrap_crosssection__array__area_cum(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%area_cum)) then
        dshape(1:1) = shape(this_ptr%p%area_cum)
        dloc = loc(this_ptr%p%area_cum)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__area_cum

subroutine f90wrap_crosssection__array__perim_cum(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%perim_cum)) then
        dshape(1:1) = shape(this_ptr%p%perim_cum)
        dloc = loc(this_ptr%p%perim_cum)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__perim_cum

subroutine f90wrap_crosssection__array__pa_cum(this, nd, dtype, dshape, dloc)
    use m_mesh, only: crosssection
    implicit none
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: this(2)
    type(crosssection_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pa_cum)) then
        dshape(1:1) = shape(this_ptr%p%pa_cum)
        dloc = loc(this_ptr%p%pa_cum)
    else
        dloc = 0
    end if
end subroutine f90wrap_crosssection__array__pa_cum

subroutine f90wrap_crosssection_initialise(cs, nlevels, shape_model)
    use m_mesh, only: crosssection, crosssection_initialise
    implicit none
    
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type(crosssection_ptr_type) :: cs_ptr
    integer, intent(out), dimension(2) :: cs
    integer(4), intent(in) :: nlevels
    character(128), intent(in), optional :: shape_model
    allocate(cs_ptr%p)
    call crosssection_initialise(cs=cs_ptr%p, nlevels=nlevels, shape_model=shape_model)
    cs = transfer(cs_ptr, cs)
end subroutine f90wrap_crosssection_initialise

subroutine f90wrap_crosssection_finalise(cs)
    use m_mesh, only: crosssection, crosssection_finalise
    implicit none
    
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type(crosssection_ptr_type) :: cs_ptr
    integer, intent(in), dimension(2) :: cs
    cs_ptr = transfer(cs, cs_ptr)
    call crosssection_finalise(cs=cs_ptr%p)
    deallocate(cs_ptr%p)
end subroutine f90wrap_crosssection_finalise

subroutine f90wrap_set_levels(cs, heights, widths, shape_model, n0, n1)
    use m_mesh, only: set_levels, crosssection
    implicit none
    
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type(crosssection_ptr_type) :: cs_ptr
    integer, intent(in), dimension(2) :: cs
    real(8), intent(in), dimension(n0) :: heights
    real(8), intent(in), dimension(n1) :: widths
    character(128), intent(in), optional :: shape_model
    integer :: n0
    !f2py intent(hide), depend(heights) :: n0 = shape(heights,0)
    integer :: n1
    !f2py intent(hide), depend(widths) :: n1 = shape(widths,0)
    cs_ptr = transfer(cs, cs_ptr)
    call set_levels(cs=cs_ptr%p, heights=heights, widths=widths, shape_model=shape_model)
end subroutine f90wrap_set_levels

subroutine f90wrap_crosssection_copy(cs_src, cs_dst)
    use m_mesh, only: crosssection_copy, crosssection
    implicit none
    
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type(crosssection_ptr_type) :: cs_src_ptr
    integer, intent(in), dimension(2) :: cs_src
    type(crosssection_ptr_type) :: cs_dst_ptr
    integer, intent(in), dimension(2) :: cs_dst
    cs_src_ptr = transfer(cs_src, cs_src_ptr)
    cs_dst_ptr = transfer(cs_dst, cs_dst_ptr)
    call crosssection_copy(cs_src=cs_src_ptr%p, cs_dst=cs_dst_ptr%p)
end subroutine f90wrap_crosssection_copy

subroutine f90wrap_update_geometry(cs)
    use m_mesh, only: update_geometry, crosssection
    implicit none
    
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    type(crosssection_ptr_type) :: cs_ptr
    integer, intent(in), dimension(2) :: cs
    cs_ptr = transfer(cs, cs_ptr)
    call update_geometry(cs=cs_ptr%p)
end subroutine f90wrap_update_geometry

subroutine f90wrap_segment__get__first_cs(this, f90wrap_first_cs)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_first_cs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_first_cs = this_ptr%p%first_cs
end subroutine f90wrap_segment__get__first_cs

subroutine f90wrap_segment__set__first_cs(this, f90wrap_first_cs)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_first_cs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%first_cs = f90wrap_first_cs
end subroutine f90wrap_segment__set__first_cs

subroutine f90wrap_segment__get__last_cs(this, f90wrap_last_cs)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_last_cs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_last_cs = this_ptr%p%last_cs
end subroutine f90wrap_segment__get__last_cs

subroutine f90wrap_segment__set__last_cs(this, f90wrap_last_cs)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_last_cs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%last_cs = f90wrap_last_cs
end subroutine f90wrap_segment__set__last_cs

subroutine f90wrap_segment__get__ds_seg(this, f90wrap_ds_seg)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ds_seg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ds_seg = this_ptr%p%ds_seg
end subroutine f90wrap_segment__get__ds_seg

subroutine f90wrap_segment__set__ds_seg(this, f90wrap_ds_seg)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ds_seg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ds_seg = f90wrap_ds_seg
end subroutine f90wrap_segment__set__ds_seg

subroutine f90wrap_segment__array__us_seg(this, nd, dtype, dshape, dloc)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in) :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%us_seg)) then
        dshape(1:1) = shape(this_ptr%p%us_seg)
        dloc = loc(this_ptr%p%us_seg)
    else
        dloc = 0
    end if
end subroutine f90wrap_segment__array__us_seg

subroutine f90wrap_segment__get__ds_bc(this, f90wrap_ds_bc)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ds_bc
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ds_bc = this_ptr%p%ds_bc
end subroutine f90wrap_segment__get__ds_bc

subroutine f90wrap_segment__set__ds_bc(this, f90wrap_ds_bc)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ds_bc
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ds_bc = f90wrap_ds_bc
end subroutine f90wrap_segment__set__ds_bc

subroutine f90wrap_segment__get__us_bc(this, f90wrap_us_bc)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_us_bc
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_us_bc = this_ptr%p%us_bc
end subroutine f90wrap_segment__get__us_bc

subroutine f90wrap_segment__set__us_bc(this, f90wrap_us_bc)
    use m_mesh, only: segment
    implicit none
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in)   :: this(2)
    type(segment_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_us_bc
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%us_bc = f90wrap_us_bc
end subroutine f90wrap_segment__set__us_bc

subroutine f90wrap_segment_initialise(seg, nus_segs)
    use m_mesh, only: segment_initialise, segment
    implicit none
    
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    type(segment_ptr_type) :: seg_ptr
    integer, intent(out), dimension(2) :: seg
    integer, intent(in) :: nus_segs
    allocate(seg_ptr%p)
    call segment_initialise(seg=seg_ptr%p, nus_segs=nus_segs)
    seg = transfer(seg_ptr, seg)
end subroutine f90wrap_segment_initialise

subroutine f90wrap_segment_finalise(seg)
    use m_mesh, only: segment, segment_finalise
    implicit none
    
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    type(segment_ptr_type) :: seg_ptr
    integer, intent(in), dimension(2) :: seg
    seg_ptr = transfer(seg, seg_ptr)
    call segment_finalise(seg=seg_ptr%p)
    deallocate(seg_ptr%p)
end subroutine f90wrap_segment_finalise

subroutine f90wrap_spatialfield__get__interp(this, f90wrap_interp)
    use m_mesh, only: spatialfield
    implicit none
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in)   :: this(2)
    type(spatialfield_ptr_type) :: this_ptr
    character(8), intent(out) :: f90wrap_interp
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_interp = this_ptr%p%interp
end subroutine f90wrap_spatialfield__get__interp

subroutine f90wrap_spatialfield__set__interp(this, f90wrap_interp)
    use m_mesh, only: spatialfield
    implicit none
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in)   :: this(2)
    type(spatialfield_ptr_type) :: this_ptr
    character(8), intent(in) :: f90wrap_interp
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%interp = f90wrap_interp
end subroutine f90wrap_spatialfield__set__interp

subroutine f90wrap_spatialfield__array__x(this, nd, dtype, dshape, dloc)
    use m_mesh, only: spatialfield
    implicit none
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in) :: this(2)
    type(spatialfield_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x)) then
        dshape(1:1) = shape(this_ptr%p%x)
        dloc = loc(this_ptr%p%x)
    else
        dloc = 0
    end if
end subroutine f90wrap_spatialfield__array__x

subroutine f90wrap_spatialfield__array__y(this, nd, dtype, dshape, dloc)
    use m_mesh, only: spatialfield
    implicit none
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in) :: this(2)
    type(spatialfield_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y)) then
        dshape(1:1) = shape(this_ptr%p%y)
        dloc = loc(this_ptr%p%y)
    else
        dloc = 0
    end if
end subroutine f90wrap_spatialfield__array__y

subroutine f90wrap_spatialfield_initialise(this)
    use m_mesh, only: spatialfield
    implicit none
    
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    type(spatialfield_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_spatialfield_initialise

subroutine f90wrap_spatialfield_finalise(this)
    use m_mesh, only: spatialfield
    implicit none
    
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    type(spatialfield_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_spatialfield_finalise

subroutine f90wrap_mesh__get__ncs(this, f90wrap_ncs)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ncs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ncs = this_ptr%p%ncs
end subroutine f90wrap_mesh__get__ncs

subroutine f90wrap_mesh__set__ncs(this, f90wrap_ncs)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ncs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncs = f90wrap_ncs
end subroutine f90wrap_mesh__set__ncs

subroutine f90wrap_mesh__get__nseg(this, f90wrap_nseg)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nseg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nseg = this_ptr%p%nseg
end subroutine f90wrap_mesh__get__nseg

subroutine f90wrap_mesh__set__nseg(this, f90wrap_nseg)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nseg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nseg = f90wrap_nseg
end subroutine f90wrap_mesh__set__nseg

subroutine f90wrap_mesh__array_getitem__cs(f90wrap_this, f90wrap_i, csitem)
    
    use m_mesh, only: mesh, crosssection
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: csitem(2)
    type(crosssection_ptr_type) :: cs_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cs)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cs)) then
            call f90wrap_abort("array index out of range")
        else
            cs_ptr%p => this_ptr%p%cs(f90wrap_i)
            csitem = transfer(cs_ptr,csitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_getitem__cs

subroutine f90wrap_mesh__array_setitem__cs(f90wrap_this, f90wrap_i, csitem)
    
    use m_mesh, only: mesh, crosssection
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: csitem(2)
    type(crosssection_ptr_type) :: cs_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cs)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cs)) then
            call f90wrap_abort("array index out of range")
        else
            cs_ptr = transfer(csitem,cs_ptr)
            this_ptr%p%cs(f90wrap_i) = cs_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_setitem__cs

subroutine f90wrap_mesh__array_len__cs(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: mesh, crosssection
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type crosssection_ptr_type
        type(crosssection), pointer :: p => NULL()
    end type crosssection_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cs)) then
        f90wrap_n = size(this_ptr%p%cs)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_mesh__array_len__cs

subroutine f90wrap_mesh__array_getitem__seg(f90wrap_this, f90wrap_i, segitem)
    
    use m_mesh, only: mesh, segment
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: segitem(2)
    type(segment_ptr_type) :: seg_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%seg)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%seg)) then
            call f90wrap_abort("array index out of range")
        else
            seg_ptr%p => this_ptr%p%seg(f90wrap_i)
            segitem = transfer(seg_ptr,segitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_getitem__seg

subroutine f90wrap_mesh__array_setitem__seg(f90wrap_this, f90wrap_i, segitem)
    
    use m_mesh, only: mesh, segment
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: segitem(2)
    type(segment_ptr_type) :: seg_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%seg)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%seg)) then
            call f90wrap_abort("array index out of range")
        else
            seg_ptr = transfer(segitem,seg_ptr)
            this_ptr%p%seg(f90wrap_i) = seg_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_setitem__seg

subroutine f90wrap_mesh__array_len__seg(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: mesh, segment
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type segment_ptr_type
        type(segment), pointer :: p => NULL()
    end type segment_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%seg)) then
        f90wrap_n = size(this_ptr%p%seg)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_mesh__array_len__seg

subroutine f90wrap_mesh__get__bathy_field(this, f90wrap_bathy_field)
    use m_mesh, only: mesh, spatialfield
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_bathy_field(2)
    type(spatialfield_ptr_type) :: bathy_field_ptr
    
    this_ptr = transfer(this, this_ptr)
    bathy_field_ptr%p => this_ptr%p%bathy_field
    f90wrap_bathy_field = transfer(bathy_field_ptr,f90wrap_bathy_field)
end subroutine f90wrap_mesh__get__bathy_field

subroutine f90wrap_mesh__set__bathy_field(this, f90wrap_bathy_field)
    use m_mesh, only: mesh, spatialfield
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_bathy_field(2)
    type(spatialfield_ptr_type) :: bathy_field_ptr
    
    this_ptr = transfer(this, this_ptr)
    bathy_field_ptr = transfer(f90wrap_bathy_field,bathy_field_ptr)
    this_ptr%p%bathy_field = bathy_field_ptr%p
end subroutine f90wrap_mesh__set__bathy_field

subroutine f90wrap_mesh__get__strickler_type(this, f90wrap_strickler_type)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    character(16), intent(out) :: f90wrap_strickler_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_strickler_type = this_ptr%p%strickler_type
end subroutine f90wrap_mesh__get__strickler_type

subroutine f90wrap_mesh__set__strickler_type(this, f90wrap_strickler_type)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    character(16), intent(in) :: f90wrap_strickler_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%strickler_type = f90wrap_strickler_type
end subroutine f90wrap_mesh__set__strickler_type

subroutine f90wrap_mesh__get__strickler_type_code(this, f90wrap_strickler_type_code)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_strickler_type_code
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_strickler_type_code = this_ptr%p%strickler_type_code
end subroutine f90wrap_mesh__get__strickler_type_code

subroutine f90wrap_mesh__set__strickler_type_code(this, f90wrap_strickler_type_code)
    use m_mesh, only: mesh
    implicit none
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_strickler_type_code
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%strickler_type_code = f90wrap_strickler_type_code
end subroutine f90wrap_mesh__set__strickler_type_code

subroutine f90wrap_mesh__array_getitem__strickler_fields(f90wrap_this, f90wrap_i, strickler_fieldsitem)
    
    use m_mesh, only: mesh, spatialfield
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: strickler_fieldsitem(2)
    type(spatialfield_ptr_type) :: strickler_fields_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%strickler_fields)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%strickler_fields)) then
            call f90wrap_abort("array index out of range")
        else
            strickler_fields_ptr%p => this_ptr%p%strickler_fields(f90wrap_i)
            strickler_fieldsitem = transfer(strickler_fields_ptr,strickler_fieldsitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_getitem__strickler_fields

subroutine f90wrap_mesh__array_setitem__strickler_fields(f90wrap_this, f90wrap_i, strickler_fieldsitem)
    
    use m_mesh, only: mesh, spatialfield
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: strickler_fieldsitem(2)
    type(spatialfield_ptr_type) :: strickler_fields_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%strickler_fields)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%strickler_fields)) then
            call f90wrap_abort("array index out of range")
        else
            strickler_fields_ptr = transfer(strickler_fieldsitem,strickler_fields_ptr)
            this_ptr%p%strickler_fields(f90wrap_i) = strickler_fields_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_mesh__array_setitem__strickler_fields

subroutine f90wrap_mesh__array_len__strickler_fields(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: mesh, spatialfield
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type spatialfield_ptr_type
        type(spatialfield), pointer :: p => NULL()
    end type spatialfield_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(mesh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%strickler_fields)) then
        f90wrap_n = size(this_ptr%p%strickler_fields)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_mesh__array_len__strickler_fields

subroutine f90wrap_dealloc_mesh(msh)
    use m_mesh, only: mesh, dealloc_mesh
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    msh_ptr = transfer(msh, msh_ptr)
    call dealloc_mesh(msh=msh_ptr%p)
end subroutine f90wrap_dealloc_mesh

subroutine f90wrap_mesh_initialise(msh, ncs, nseg)
    use m_mesh, only: mesh, mesh_initialise
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(out), dimension(2) :: msh
    integer, intent(in) :: ncs
    integer, intent(in), optional :: nseg
    allocate(msh_ptr%p)
    call mesh_initialise(msh=msh_ptr%p, ncs=ncs, nseg=nseg)
    msh = transfer(msh_ptr, msh)
end subroutine f90wrap_mesh_initialise

subroutine f90wrap_setup_segment(msh, iseg, first_cs, ncs, us_seg, ds_seg, n0)
    use m_mesh, only: mesh, setup_segment
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    integer(4), intent(in) :: iseg
    integer(4), intent(in) :: first_cs
    integer(4), intent(in) :: ncs
    integer(4), intent(in), dimension(n0) :: us_seg
    integer(4), intent(in) :: ds_seg
    integer :: n0
    !f2py intent(hide), depend(us_seg) :: n0 = shape(us_seg,0)
    msh_ptr = transfer(msh, msh_ptr)
    call setup_segment(msh=msh_ptr%p, iseg=iseg, first_cs=first_cs, ncs=ncs, us_seg=us_seg, ds_seg=ds_seg)
end subroutine f90wrap_setup_segment

subroutine f90wrap_setup_crosssection(msh, ics, nlevels)
    use m_mesh, only: mesh, setup_crosssection
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    integer(4), intent(in) :: ics
    integer(4), intent(in) :: nlevels
    msh_ptr = transfer(msh, msh_ptr)
    call setup_crosssection(msh=msh_ptr%p, ics=ics, nlevels=nlevels)
end subroutine f90wrap_setup_crosssection

subroutine f90wrap_update_geometries(msh)
    use m_mesh, only: mesh, update_geometries
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    msh_ptr = transfer(msh, msh_ptr)
    call update_geometries(msh=msh_ptr%p)
end subroutine f90wrap_update_geometries

subroutine f90wrap_set_strickler_type(msh, strickler_type)
    use m_mesh, only: mesh, set_strickler_type
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    character(*), intent(in) :: strickler_type
    msh_ptr = transfer(msh, msh_ptr)
    call set_strickler_type(msh=msh_ptr%p, strickler_type=strickler_type)
end subroutine f90wrap_set_strickler_type

subroutine f90wrap_set_uniform_strickler_parameters(msh, strickler_params)
    use m_mesh, only: mesh, set_uniform_strickler_parameters
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    real(8), dimension(2), intent(in) :: strickler_params
    msh_ptr = transfer(msh, msh_ptr)
    call set_uniform_strickler_parameters(msh=msh_ptr%p, strickler_params=strickler_params)
end subroutine f90wrap_set_uniform_strickler_parameters

subroutine f90wrap_set_bathy_field_linear(msh, x, bathy, n0, n1)
    use m_mesh, only: mesh, set_bathy_field_linear
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(in), dimension(n1) :: bathy
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(bathy) :: n1 = shape(bathy,0)
    msh_ptr = transfer(msh, msh_ptr)
    call set_bathy_field_linear(msh=msh_ptr%p, x=x, bathy=bathy)
end subroutine f90wrap_set_bathy_field_linear

subroutine f90wrap_set_strickler_fields_segment(msh, strickler_params, n0, n1)
    use m_mesh, only: set_strickler_fields_segment, mesh
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    real(8), intent(in), dimension(n0,n1) :: strickler_params
    integer :: n0
    !f2py intent(hide), depend(strickler_params) :: n0 = shape(strickler_params,0)
    integer :: n1
    !f2py intent(hide), depend(strickler_params) :: n1 = shape(strickler_params,1)
    msh_ptr = transfer(msh, msh_ptr)
    call set_strickler_fields_segment(msh=msh_ptr%p, strickler_params=strickler_params)
end subroutine f90wrap_set_strickler_fields_segment

subroutine f90wrap_set_strickler_fields_linear(msh, x, strickler_params, n0, n1, n2)
    use m_mesh, only: mesh, set_strickler_fields_linear
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(in), dimension(n1,n2) :: strickler_params
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(strickler_params) :: n1 = shape(strickler_params,0)
    integer :: n2
    !f2py intent(hide), depend(strickler_params) :: n2 = shape(strickler_params,1)
    msh_ptr = transfer(msh, msh_ptr)
    call set_strickler_fields_linear(msh=msh_ptr%p, x=x, strickler_params=strickler_params)
end subroutine f90wrap_set_strickler_fields_linear

subroutine f90wrap_get_segment_cs_index(msh, iseg, ics_seg, ics)
    use m_mesh, only: mesh, get_segment_cs_index
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    integer(4), intent(in) :: iseg
    integer(4), intent(in) :: ics_seg
    integer(4), intent(out) :: ics
    msh_ptr = transfer(msh, msh_ptr)
    call get_segment_cs_index(msh=msh_ptr%p, iseg=iseg, ics_seg=ics_seg, ics=ics)
end subroutine f90wrap_get_segment_cs_index

subroutine f90wrap_mesh_finalise(this)
    use m_mesh, only: mesh
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(mesh_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_mesh_finalise

subroutine f90wrap_point_in_mesh__get__indexi(this, f90wrap_indexi)
    use m_mesh, only: point_in_mesh
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_indexi
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_indexi = this_ptr%p%indexi
end subroutine f90wrap_point_in_mesh__get__indexi

subroutine f90wrap_point_in_mesh__set__indexi(this, f90wrap_indexi)
    use m_mesh, only: point_in_mesh
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_indexi
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%indexi = f90wrap_indexi
end subroutine f90wrap_point_in_mesh__set__indexi

subroutine f90wrap_point_in_mesh_initialise(this)
    use m_mesh, only: point_in_mesh
    implicit none
    
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_point_in_mesh_initialise

subroutine f90wrap_point_in_mesh_finalise(this)
    use m_mesh, only: point_in_mesh
    implicit none
    
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_point_in_mesh_finalise

subroutine f90wrap_m_mesh__get__K_constant(f90wrap_K_constant)
    use m_mesh, only: m_mesh_K_constant => K_constant
    implicit none
    integer, intent(out) :: f90wrap_K_constant
    
    f90wrap_K_constant = m_mesh_K_constant
end subroutine f90wrap_m_mesh__get__K_constant

subroutine f90wrap_m_mesh__get__K_powerlaw_h(f90wrap_K_powerlaw_h)
    use m_mesh, only: m_mesh_K_powerlaw_h => K_powerlaw_h
    implicit none
    integer, intent(out) :: f90wrap_K_powerlaw_h
    
    f90wrap_K_powerlaw_h = m_mesh_K_powerlaw_h
end subroutine f90wrap_m_mesh__get__K_powerlaw_h

subroutine f90wrap_m_mesh__get__strickler_type_constant(f90wrap_strickler_type_constant)
    use m_mesh, only: m_mesh_strickler_type_constant => strickler_type_constant
    implicit none
    integer, intent(out) :: f90wrap_strickler_type_constant
    
    f90wrap_strickler_type_constant = m_mesh_strickler_type_constant
end subroutine f90wrap_m_mesh__get__strickler_type_constant

subroutine f90wrap_m_mesh__get__strickler_type_powerlaw_h(f90wrap_strickler_type_powerlaw_h)
    use m_mesh, only: m_mesh_strickler_type_powerlaw_h => strickler_type_powerlaw_h
    implicit none
    integer, intent(out) :: f90wrap_strickler_type_powerlaw_h
    
    f90wrap_strickler_type_powerlaw_h = m_mesh_strickler_type_powerlaw_h
end subroutine f90wrap_m_mesh__get__strickler_type_powerlaw_h

! End of module m_mesh defined in file m_mesh.f90

