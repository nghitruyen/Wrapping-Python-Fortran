! Module m_sw_mono defined in file m_sw_mono.f90

subroutine f90wrap_unknowns__array__A(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%A)) then
        dshape(1:1) = shape(this_ptr%p%A)
        dloc = loc(this_ptr%p%A)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__A

subroutine f90wrap_unknowns__array__Q(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Q)) then
        dshape(1:1) = shape(this_ptr%p%Q)
        dloc = loc(this_ptr%p%Q)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__Q

subroutine f90wrap_unknowns__array__h(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__h

subroutine f90wrap_unknowns__array__qlat(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%qlat)) then
        dshape(1:2) = shape(this_ptr%p%qlat)
        dloc = loc(this_ptr%p%qlat)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__qlat

subroutine f90wrap_unknowns__array__Sf(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Sf)) then
        dshape(1:1) = shape(this_ptr%p%Sf)
        dloc = loc(this_ptr%p%Sf)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__Sf

subroutine f90wrap_unknowns__array__Sg(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: unknowns
    implicit none
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in) :: this(2)
    type(unknowns_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%Sg)) then
        dshape(1:1) = shape(this_ptr%p%Sg)
        dloc = loc(this_ptr%p%Sg)
    else
        dloc = 0
    end if
end subroutine f90wrap_unknowns__array__Sg

subroutine f90wrap_unknowns_initialise(dof, msh)
    use m_mesh, only: mesh
    use m_sw_mono, only: unknowns_initialise, unknowns
    implicit none
    
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    type(unknowns_ptr_type) :: dof_ptr
    integer, intent(out), dimension(2) :: dof
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    msh_ptr = transfer(msh, msh_ptr)
    allocate(dof_ptr%p)
    call unknowns_initialise(dof=dof_ptr%p, msh=msh_ptr%p)
    dof = transfer(dof_ptr, dof)
end subroutine f90wrap_unknowns_initialise

subroutine f90wrap_unknowns_finalise(dof)
    use m_sw_mono, only: unknowns_finalise, unknowns
    implicit none
    
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    type(unknowns_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    dof_ptr = transfer(dof, dof_ptr)
    call unknowns_finalise(dof=dof_ptr%p)
    deallocate(dof_ptr%p)
end subroutine f90wrap_unknowns_finalise

subroutine f90wrap_implicitmatrix__array__seg_offsets(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%seg_offsets)) then
        dshape(1:1) = shape(this_ptr%p%seg_offsets)
        dloc = loc(this_ptr%p%seg_offsets)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__seg_offsets

subroutine f90wrap_implicitmatrix__array__GA(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GA)) then
        dshape(1:1) = shape(this_ptr%p%GA)
        dloc = loc(this_ptr%p%GA)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GA

subroutine f90wrap_implicitmatrix__array__GB(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GB)) then
        dshape(1:1) = shape(this_ptr%p%GB)
        dloc = loc(this_ptr%p%GB)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GB

subroutine f90wrap_implicitmatrix__array__GC(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GC)) then
        dshape(1:1) = shape(this_ptr%p%GC)
        dloc = loc(this_ptr%p%GC)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GC

subroutine f90wrap_implicitmatrix__array__GD(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GD)) then
        dshape(1:1) = shape(this_ptr%p%GD)
        dloc = loc(this_ptr%p%GD)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GD

subroutine f90wrap_implicitmatrix__array__GE(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GE)) then
        dshape(1:1) = shape(this_ptr%p%GE)
        dloc = loc(this_ptr%p%GE)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GE

subroutine f90wrap_implicitmatrix__array__GF(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%GF)) then
        dshape(1:1) = shape(this_ptr%p%GF)
        dloc = loc(this_ptr%p%GF)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__GF

subroutine f90wrap_implicitmatrix__array__CR(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%CR)) then
        dshape(1:1) = shape(this_ptr%p%CR)
        dloc = loc(this_ptr%p%CR)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__CR

subroutine f90wrap_implicitmatrix__array__CS(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%CS)) then
        dshape(1:1) = shape(this_ptr%p%CS)
        dloc = loc(this_ptr%p%CS)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__CS

subroutine f90wrap_implicitmatrix__array__CT(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%CT)) then
        dshape(1:1) = shape(this_ptr%p%CT)
        dloc = loc(this_ptr%p%CT)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__CT

subroutine f90wrap_implicitmatrix__array__TA1(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TA1)) then
        dshape(1:1) = shape(this_ptr%p%TA1)
        dloc = loc(this_ptr%p%TA1)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TA1

subroutine f90wrap_implicitmatrix__array__TA2(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TA2)) then
        dshape(1:1) = shape(this_ptr%p%TA2)
        dloc = loc(this_ptr%p%TA2)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TA2

subroutine f90wrap_implicitmatrix__array__TA3(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TA3)) then
        dshape(1:1) = shape(this_ptr%p%TA3)
        dloc = loc(this_ptr%p%TA3)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TA3

subroutine f90wrap_implicitmatrix__array__TA4(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TA4)) then
        dshape(1:1) = shape(this_ptr%p%TA4)
        dloc = loc(this_ptr%p%TA4)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TA4

subroutine f90wrap_implicitmatrix__array__TB1(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TB1)) then
        dshape(1:1) = shape(this_ptr%p%TB1)
        dloc = loc(this_ptr%p%TB1)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TB1

subroutine f90wrap_implicitmatrix__array__TB2(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%TB2)) then
        dshape(1:1) = shape(this_ptr%p%TB2)
        dloc = loc(this_ptr%p%TB2)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__TB2

subroutine f90wrap_implicitmatrix__array__ANZ(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ANZ)) then
        dshape(1:1) = shape(this_ptr%p%ANZ)
        dloc = loc(this_ptr%p%ANZ)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__ANZ

subroutine f90wrap_implicitmatrix__array__RHS(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: implicitmatrix
    implicit none
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in) :: this(2)
    type(implicitmatrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%RHS)) then
        dshape(1:1) = shape(this_ptr%p%RHS)
        dloc = loc(this_ptr%p%RHS)
    else
        dloc = 0
    end if
end subroutine f90wrap_implicitmatrix__array__RHS

subroutine f90wrap_implicitmatrix_initialise(imp, mdl, msh)
    use m_sw_mono, only: model, implicitmatrix, implicitmatrix_initialise
    use m_mesh, only: mesh
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    type(implicitmatrix_ptr_type) :: imp_ptr
    integer, intent(out), dimension(2) :: imp
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    mdl_ptr = transfer(mdl, mdl_ptr)
    msh_ptr = transfer(msh, msh_ptr)
    allocate(imp_ptr%p)
    call implicitmatrix_initialise(imp=imp_ptr%p, mdl=mdl_ptr%p, msh=msh_ptr%p)
    imp = transfer(imp_ptr, imp)
end subroutine f90wrap_implicitmatrix_initialise

subroutine f90wrap_implicitmatrix_finalise(imp)
    use m_sw_mono, only: implicitmatrix_finalise, implicitmatrix
    implicit none
    
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    type(implicitmatrix_ptr_type) :: imp_ptr
    integer, intent(in), dimension(2) :: imp
    imp_ptr = transfer(imp, imp_ptr)
    call implicitmatrix_finalise(imp=imp_ptr%p)
    deallocate(imp_ptr%p)
end subroutine f90wrap_implicitmatrix_finalise

subroutine f90wrap_timeseries__array__t(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: timeseries
    implicit none
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in) :: this(2)
    type(timeseries_ptr_type) :: this_ptr
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
end subroutine f90wrap_timeseries__array__t

subroutine f90wrap_timeseries__array__y(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: timeseries
    implicit none
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in) :: this(2)
    type(timeseries_ptr_type) :: this_ptr
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
end subroutine f90wrap_timeseries__array__y

subroutine f90wrap_timeseries_initialise(ts, nt, y0)
    use m_sw_mono, only: timeseries_initialise, timeseries
    implicit none
    
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    type(timeseries_ptr_type) :: ts_ptr
    integer, intent(out), dimension(2) :: ts
    integer(4), intent(in) :: nt
    real(8), optional :: y0
    allocate(ts_ptr%p)
    call timeseries_initialise(ts=ts_ptr%p, nt=nt, y0=y0)
    ts = transfer(ts_ptr, ts)
end subroutine f90wrap_timeseries_initialise

subroutine f90wrap_timeseries_finalise(ts)
    use m_sw_mono, only: timeseries, timeseries_finalise
    implicit none
    
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    type(timeseries_ptr_type) :: ts_ptr
    integer, intent(in), dimension(2) :: ts
    ts_ptr = transfer(ts, ts_ptr)
    call timeseries_finalise(ts=ts_ptr%p)
    deallocate(ts_ptr%p)
end subroutine f90wrap_timeseries_finalise

subroutine f90wrap_boundarycondition__get__id(this, f90wrap_id)
    use m_sw_mono, only: boundarycondition
    implicit none
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    integer, intent(in)   :: this(2)
    type(boundarycondition_ptr_type) :: this_ptr
    character(16), intent(out) :: f90wrap_id
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_id = this_ptr%p%id
end subroutine f90wrap_boundarycondition__get__id

subroutine f90wrap_boundarycondition__set__id(this, f90wrap_id)
    use m_sw_mono, only: boundarycondition
    implicit none
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    integer, intent(in)   :: this(2)
    type(boundarycondition_ptr_type) :: this_ptr
    character(16), intent(in) :: f90wrap_id
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%id = f90wrap_id
end subroutine f90wrap_boundarycondition__set__id

subroutine f90wrap_boundarycondition__get__ts(this, f90wrap_ts)
    use m_sw_mono, only: timeseries, boundarycondition
    implicit none
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in)   :: this(2)
    type(boundarycondition_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ts(2)
    type(timeseries_ptr_type) :: ts_ptr
    
    this_ptr = transfer(this, this_ptr)
    ts_ptr%p => this_ptr%p%ts
    f90wrap_ts = transfer(ts_ptr,f90wrap_ts)
end subroutine f90wrap_boundarycondition__get__ts

subroutine f90wrap_boundarycondition__set__ts(this, f90wrap_ts)
    use m_sw_mono, only: timeseries, boundarycondition
    implicit none
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in)   :: this(2)
    type(boundarycondition_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ts(2)
    type(timeseries_ptr_type) :: ts_ptr
    
    this_ptr = transfer(this, this_ptr)
    ts_ptr = transfer(f90wrap_ts,ts_ptr)
    this_ptr%p%ts = ts_ptr%p
end subroutine f90wrap_boundarycondition__set__ts

subroutine f90wrap_set_timeseries(bc, t, y, n0, n1)
    use m_sw_mono, only: boundarycondition, set_timeseries
    implicit none
    
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    type(boundarycondition_ptr_type) :: bc_ptr
    integer, intent(in), dimension(2) :: bc
    real(8), intent(in), dimension(n0) :: t
    real(8), intent(in), dimension(n1) :: y
    integer :: n0
    !f2py intent(hide), depend(t) :: n0 = shape(t,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    bc_ptr = transfer(bc, bc_ptr)
    call set_timeseries(bc=bc_ptr%p, t=t, y=y)
end subroutine f90wrap_set_timeseries

subroutine f90wrap_boundarycondition_initialise(this)
    use m_sw_mono, only: boundarycondition
    implicit none
    
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    type(boundarycondition_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_boundarycondition_initialise

subroutine f90wrap_boundarycondition_finalise(this)
    use m_sw_mono, only: boundarycondition
    implicit none
    
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    type(boundarycondition_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_boundarycondition_finalise

subroutine f90wrap_inflowcondition__get__iseg(this, f90wrap_iseg)
    use m_sw_mono, only: inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_iseg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_iseg = this_ptr%p%iseg
end subroutine f90wrap_inflowcondition__get__iseg

subroutine f90wrap_inflowcondition__set__iseg(this, f90wrap_iseg)
    use m_sw_mono, only: inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_iseg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%iseg = f90wrap_iseg
end subroutine f90wrap_inflowcondition__set__iseg

subroutine f90wrap_inflowcondition__get__ie(this, f90wrap_ie)
    use m_sw_mono, only: inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ie
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ie = this_ptr%p%ie
end subroutine f90wrap_inflowcondition__get__ie

subroutine f90wrap_inflowcondition__set__ie(this, f90wrap_ie)
    use m_sw_mono, only: inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ie
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ie = f90wrap_ie
end subroutine f90wrap_inflowcondition__set__ie

subroutine f90wrap_inflowcondition__get__ts(this, f90wrap_ts)
    use m_sw_mono, only: timeseries, inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_ts(2)
    type(timeseries_ptr_type) :: ts_ptr
    
    this_ptr = transfer(this, this_ptr)
    ts_ptr%p => this_ptr%p%ts
    f90wrap_ts = transfer(ts_ptr,f90wrap_ts)
end subroutine f90wrap_inflowcondition__get__ts

subroutine f90wrap_inflowcondition__set__ts(this, f90wrap_ts)
    use m_sw_mono, only: timeseries, inflowcondition
    implicit none
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    type timeseries_ptr_type
        type(timeseries), pointer :: p => NULL()
    end type timeseries_ptr_type
    integer, intent(in)   :: this(2)
    type(inflowcondition_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_ts(2)
    type(timeseries_ptr_type) :: ts_ptr
    
    this_ptr = transfer(this, this_ptr)
    ts_ptr = transfer(f90wrap_ts,ts_ptr)
    this_ptr%p%ts = ts_ptr%p
end subroutine f90wrap_inflowcondition__set__ts

subroutine f90wrap_inflowcondition_initialise(this)
    use m_sw_mono, only: inflowcondition
    implicit none
    
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    type(inflowcondition_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_inflowcondition_initialise

subroutine f90wrap_inflowcondition_finalise(this)
    use m_sw_mono, only: inflowcondition
    implicit none
    
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    type(inflowcondition_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_inflowcondition_finalise

subroutine f90wrap_results__array__t(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: results
    implicit none
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in) :: this(2)
    type(results_ptr_type) :: this_ptr
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
end subroutine f90wrap_results__array__t

subroutine f90wrap_results__array__q(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: results
    implicit none
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in) :: this(2)
    type(results_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:2) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_results__array__q

subroutine f90wrap_results__array__h(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: results
    implicit none
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in) :: this(2)
    type(results_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:2) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_results__array__h

subroutine f90wrap_results__array__a(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: results
    implicit none
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in) :: this(2)
    type(results_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%a)) then
        dshape(1:2) = shape(this_ptr%p%a)
        dloc = loc(this_ptr%p%a)
    else
        dloc = 0
    end if
end subroutine f90wrap_results__array__a

subroutine f90wrap_results_initialise(this)
    use m_sw_mono, only: results
    implicit none
    
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    type(results_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_results_initialise

subroutine f90wrap_results_finalise(this)
    use m_sw_mono, only: results
    implicit none
    
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    type(results_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_results_finalise

subroutine f90wrap_model__get__discharge_estimation(this, f90wrap_discharge_estimation)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_discharge_estimation
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_discharge_estimation = this_ptr%p%discharge_estimation
end subroutine f90wrap_model__get__discharge_estimation

subroutine f90wrap_model__set__discharge_estimation(this, f90wrap_discharge_estimation)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_discharge_estimation
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%discharge_estimation = f90wrap_discharge_estimation
end subroutine f90wrap_model__set__discharge_estimation

subroutine f90wrap_model__get__status(this, f90wrap_status)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_status
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_status = this_ptr%p%status
end subroutine f90wrap_model__get__status

subroutine f90wrap_model__set__status(this, f90wrap_status)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_status
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%status = f90wrap_status
end subroutine f90wrap_model__set__status

subroutine f90wrap_model__array__warning_counters(this, nd, dtype, dshape, dloc)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in) :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%warning_counters)
    dloc = loc(this_ptr%p%warning_counters)
end subroutine f90wrap_model__array__warning_counters

subroutine f90wrap_model__get__tc(this, f90wrap_tc)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_tc
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_tc = this_ptr%p%tc
end subroutine f90wrap_model__get__tc

subroutine f90wrap_model__set__tc(this, f90wrap_tc)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_tc
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%tc = f90wrap_tc
end subroutine f90wrap_model__set__tc

subroutine f90wrap_model__get__te(this, f90wrap_te)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_te
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_te = this_ptr%p%te
end subroutine f90wrap_model__get__te

subroutine f90wrap_model__set__te(this, f90wrap_te)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_te
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%te = f90wrap_te
end subroutine f90wrap_model__set__te

subroutine f90wrap_model__get__ts(this, f90wrap_ts)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ts
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ts = this_ptr%p%ts
end subroutine f90wrap_model__get__ts

subroutine f90wrap_model__set__ts(this, f90wrap_ts)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ts
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ts = f90wrap_ts
end subroutine f90wrap_model__set__ts

subroutine f90wrap_model__get__dt(this, f90wrap_dt)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dt = this_ptr%p%dt
end subroutine f90wrap_model__get__dt

subroutine f90wrap_model__set__dt(this, f90wrap_dt)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dt = f90wrap_dt
end subroutine f90wrap_model__set__dt

subroutine f90wrap_model__get__dtout(this, f90wrap_dtout)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dtout
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dtout = this_ptr%p%dtout
end subroutine f90wrap_model__get__dtout

subroutine f90wrap_model__set__dtout(this, f90wrap_dtout)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dtout
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dtout = f90wrap_dtout
end subroutine f90wrap_model__set__dtout

subroutine f90wrap_model__get__frLPI(this, f90wrap_frLPI)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_frLPI
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_frLPI = this_ptr%p%frLPI
end subroutine f90wrap_model__get__frLPI

subroutine f90wrap_model__set__frLPI(this, f90wrap_frLPI)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_frLPI
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%frLPI = f90wrap_frLPI
end subroutine f90wrap_model__set__frLPI

subroutine f90wrap_model__get__gravity(this, f90wrap_gravity)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_gravity
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_gravity = this_ptr%p%gravity
end subroutine f90wrap_model__get__gravity

subroutine f90wrap_model__set__gravity(this, f90wrap_gravity)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_gravity
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%gravity = f90wrap_gravity
end subroutine f90wrap_model__set__gravity

subroutine f90wrap_model__get__heps(this, f90wrap_heps)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_heps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_heps = this_ptr%p%heps
end subroutine f90wrap_model__get__heps

subroutine f90wrap_model__set__heps(this, f90wrap_heps)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_heps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%heps = f90wrap_heps
end subroutine f90wrap_model__set__heps

subroutine f90wrap_model__get__qeps(this, f90wrap_qeps)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_qeps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qeps = this_ptr%p%qeps
end subroutine f90wrap_model__get__qeps

subroutine f90wrap_model__set__qeps(this, f90wrap_qeps)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_qeps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qeps = f90wrap_qeps
end subroutine f90wrap_model__set__qeps

subroutine f90wrap_model__get__mLPI(this, f90wrap_mLPI)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_mLPI
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mLPI = this_ptr%p%mLPI
end subroutine f90wrap_model__get__mLPI

subroutine f90wrap_model__set__mLPI(this, f90wrap_mLPI)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_mLPI
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mLPI = f90wrap_mLPI
end subroutine f90wrap_model__set__mLPI

subroutine f90wrap_model__get__theta_preissmann(this, f90wrap_theta_preissmann)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_theta_preissmann
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_theta_preissmann = this_ptr%p%theta_preissmann
end subroutine f90wrap_model__get__theta_preissmann

subroutine f90wrap_model__set__theta_preissmann(this, f90wrap_theta_preissmann)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_theta_preissmann
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%theta_preissmann = f90wrap_theta_preissmann
end subroutine f90wrap_model__set__theta_preissmann

subroutine f90wrap_model__get__gamma_reg(this, f90wrap_gamma_reg)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_gamma_reg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_gamma_reg = this_ptr%p%gamma_reg
end subroutine f90wrap_model__get__gamma_reg

subroutine f90wrap_model__set__gamma_reg(this, f90wrap_gamma_reg)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_gamma_reg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%gamma_reg = f90wrap_gamma_reg
end subroutine f90wrap_model__set__gamma_reg

subroutine f90wrap_model__get__cost_obs(this, f90wrap_cost_obs)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_cost_obs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cost_obs = this_ptr%p%cost_obs
end subroutine f90wrap_model__get__cost_obs

subroutine f90wrap_model__set__cost_obs(this, f90wrap_cost_obs)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_cost_obs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cost_obs = f90wrap_cost_obs
end subroutine f90wrap_model__set__cost_obs

subroutine f90wrap_model__get__cost_reg(this, f90wrap_cost_reg)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_cost_reg
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cost_reg = this_ptr%p%cost_reg
end subroutine f90wrap_model__get__cost_reg

subroutine f90wrap_model__set__cost_reg(this, f90wrap_cost_reg)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_cost_reg
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cost_reg = f90wrap_cost_reg
end subroutine f90wrap_model__set__cost_reg

subroutine f90wrap_model__get__scheme(this, f90wrap_scheme)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    character(32), intent(out) :: f90wrap_scheme
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_scheme = this_ptr%p%scheme
end subroutine f90wrap_model__get__scheme

subroutine f90wrap_model__set__scheme(this, f90wrap_scheme)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    character(32), intent(in) :: f90wrap_scheme
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%scheme = f90wrap_scheme
end subroutine f90wrap_model__set__scheme

subroutine f90wrap_model__get__scale_idw(this, f90wrap_scale_idw)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    character(16), intent(out) :: f90wrap_scale_idw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_scale_idw = this_ptr%p%scale_idw
end subroutine f90wrap_model__get__scale_idw

subroutine f90wrap_model__set__scale_idw(this, f90wrap_scale_idw)
    use m_sw_mono, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    character(16), intent(in) :: f90wrap_scale_idw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%scale_idw = f90wrap_scale_idw
end subroutine f90wrap_model__set__scale_idw

subroutine f90wrap_model__array_getitem__bc(f90wrap_this, f90wrap_i, bcitem)
    
    use m_sw_mono, only: model, boundarycondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: bcitem(2)
    type(boundarycondition_ptr_type) :: bc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%bc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%bc)) then
            call f90wrap_abort("array index out of range")
        else
            bc_ptr%p => this_ptr%p%bc(f90wrap_i)
            bcitem = transfer(bc_ptr,bcitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_model__array_getitem__bc

subroutine f90wrap_model__array_setitem__bc(f90wrap_this, f90wrap_i, bcitem)
    
    use m_sw_mono, only: model, boundarycondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: bcitem(2)
    type(boundarycondition_ptr_type) :: bc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%bc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%bc)) then
            call f90wrap_abort("array index out of range")
        else
            bc_ptr = transfer(bcitem,bc_ptr)
            this_ptr%p%bc(f90wrap_i) = bc_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_model__array_setitem__bc

subroutine f90wrap_model__array_len__bc(f90wrap_this, f90wrap_n)
    
    use m_sw_mono, only: model, boundarycondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type boundarycondition_ptr_type
        type(boundarycondition), pointer :: p => NULL()
    end type boundarycondition_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%bc)) then
        f90wrap_n = size(this_ptr%p%bc)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_model__array_len__bc

subroutine f90wrap_model__array_getitem__ic(f90wrap_this, f90wrap_i, icitem)
    
    use m_sw_mono, only: model, inflowcondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: icitem(2)
    type(inflowcondition_ptr_type) :: ic_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ic)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%ic)) then
            call f90wrap_abort("array index out of range")
        else
            ic_ptr%p => this_ptr%p%ic(f90wrap_i)
            icitem = transfer(ic_ptr,icitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_model__array_getitem__ic

subroutine f90wrap_model__array_setitem__ic(f90wrap_this, f90wrap_i, icitem)
    
    use m_sw_mono, only: model, inflowcondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: icitem(2)
    type(inflowcondition_ptr_type) :: ic_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ic)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%ic)) then
            call f90wrap_abort("array index out of range")
        else
            ic_ptr = transfer(icitem,ic_ptr)
            this_ptr%p%ic(f90wrap_i) = ic_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_model__array_setitem__ic

subroutine f90wrap_model__array_len__ic(f90wrap_this, f90wrap_n)
    
    use m_sw_mono, only: model, inflowcondition
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type inflowcondition_ptr_type
        type(inflowcondition), pointer :: p => NULL()
    end type inflowcondition_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(model_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ic)) then
        f90wrap_n = size(this_ptr%p%ic)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_model__array_len__ic

subroutine f90wrap_model__get__msh(this, f90wrap_msh)
    use m_sw_mono, only: model
    use m_mesh, only: mesh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_msh(2)
    type(mesh_ptr_type) :: msh_ptr
    
    this_ptr = transfer(this, this_ptr)
    msh_ptr%p => this_ptr%p%msh
    f90wrap_msh = transfer(msh_ptr,f90wrap_msh)
end subroutine f90wrap_model__get__msh

subroutine f90wrap_model__set__msh(this, f90wrap_msh)
    use m_sw_mono, only: model
    use m_mesh, only: mesh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_msh(2)
    type(mesh_ptr_type) :: msh_ptr
    
    this_ptr = transfer(this, this_ptr)
    msh_ptr = transfer(f90wrap_msh,msh_ptr)
    this_ptr%p%msh = msh_ptr%p
end subroutine f90wrap_model__set__msh

subroutine f90wrap_model__get__large_grid(this, f90wrap_large_grid)
    use m_sw_mono, only: model
    use m_mesh, only: mesh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_large_grid(2)
    type(mesh_ptr_type) :: large_grid_ptr
    
    this_ptr = transfer(this, this_ptr)
    large_grid_ptr%p => this_ptr%p%large_grid
    f90wrap_large_grid = transfer(large_grid_ptr,f90wrap_large_grid)
end subroutine f90wrap_model__get__large_grid

subroutine f90wrap_model__set__large_grid(this, f90wrap_large_grid)
    use m_sw_mono, only: model
    use m_mesh, only: mesh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_large_grid(2)
    type(mesh_ptr_type) :: large_grid_ptr
    
    this_ptr = transfer(this, this_ptr)
    large_grid_ptr = transfer(f90wrap_large_grid,large_grid_ptr)
    this_ptr%p%large_grid = large_grid_ptr%p
end subroutine f90wrap_model__set__large_grid

subroutine f90wrap_model__get__dof(this, f90wrap_dof)
    use m_sw_mono, only: model, unknowns
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_dof(2)
    type(unknowns_ptr_type) :: dof_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof_ptr%p => this_ptr%p%dof
    f90wrap_dof = transfer(dof_ptr,f90wrap_dof)
end subroutine f90wrap_model__get__dof

subroutine f90wrap_model__set__dof(this, f90wrap_dof)
    use m_sw_mono, only: model, unknowns
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unknowns_ptr_type
        type(unknowns), pointer :: p => NULL()
    end type unknowns_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_dof(2)
    type(unknowns_ptr_type) :: dof_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof_ptr = transfer(f90wrap_dof,dof_ptr)
    this_ptr%p%dof = dof_ptr%p
end subroutine f90wrap_model__set__dof

subroutine f90wrap_model__get__imp(this, f90wrap_imp)
    use m_sw_mono, only: model, implicitmatrix
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_imp(2)
    type(implicitmatrix_ptr_type) :: imp_ptr
    
    this_ptr = transfer(this, this_ptr)
    imp_ptr%p => this_ptr%p%imp
    f90wrap_imp = transfer(imp_ptr,f90wrap_imp)
end subroutine f90wrap_model__get__imp

subroutine f90wrap_model__set__imp(this, f90wrap_imp)
    use m_sw_mono, only: model, implicitmatrix
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type implicitmatrix_ptr_type
        type(implicitmatrix), pointer :: p => NULL()
    end type implicitmatrix_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_imp(2)
    type(implicitmatrix_ptr_type) :: imp_ptr
    
    this_ptr = transfer(this, this_ptr)
    imp_ptr = transfer(f90wrap_imp,imp_ptr)
    this_ptr%p%imp = imp_ptr%p
end subroutine f90wrap_model__set__imp

subroutine f90wrap_model__get__res(this, f90wrap_res)
    use m_sw_mono, only: model, results
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_res(2)
    type(results_ptr_type) :: res_ptr
    
    this_ptr = transfer(this, this_ptr)
    res_ptr%p => this_ptr%p%res
    f90wrap_res = transfer(res_ptr,f90wrap_res)
end subroutine f90wrap_model__get__res

subroutine f90wrap_model__set__res(this, f90wrap_res)
    use m_sw_mono, only: model, results
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type results_ptr_type
        type(results), pointer :: p => NULL()
    end type results_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_res(2)
    type(results_ptr_type) :: res_ptr
    
    this_ptr = transfer(this, this_ptr)
    res_ptr = transfer(f90wrap_res,res_ptr)
    this_ptr%p%res = res_ptr%p
end subroutine f90wrap_model__set__res

subroutine f90wrap_model_initialise(mdl, msh, scheme)
    use m_sw_mono, only: model, model_initialise
    use m_mesh, only: mesh
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(out), dimension(2) :: mdl
    type(mesh_ptr_type) :: msh_ptr
    integer, intent(in), dimension(2) :: msh
    character(*), intent(in), optional :: scheme
    msh_ptr = transfer(msh, msh_ptr)
    allocate(mdl_ptr%p)
    call model_initialise(mdl=mdl_ptr%p, msh=msh_ptr%p, scheme=scheme)
    mdl = transfer(mdl_ptr, mdl)
end subroutine f90wrap_model_initialise

subroutine f90wrap_set_large_grid(mdl, large_grid)
    use m_sw_mono, only: set_large_grid, model
    use m_mesh, only: mesh
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type mesh_ptr_type
        type(mesh), pointer :: p => NULL()
    end type mesh_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    type(mesh_ptr_type) :: large_grid_ptr
    integer, intent(in), dimension(2) :: large_grid
    mdl_ptr = transfer(mdl, mdl_ptr)
    large_grid_ptr = transfer(large_grid, large_grid_ptr)
    call set_large_grid(mdl=mdl_ptr%p, large_grid=large_grid_ptr%p)
end subroutine f90wrap_set_large_grid

subroutine f90wrap_add_inflow_condition(mdl, iseg, coords, t, q, n0, n1)
    use m_sw_mono, only: add_inflow_condition, model
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    integer(4) :: iseg
    real(8), dimension(2) :: coords
    real(8), dimension(n0) :: t
    real(8), dimension(n1) :: q
    integer :: n0
    !f2py intent(hide), depend(t) :: n0 = shape(t,0)
    integer :: n1
    !f2py intent(hide), depend(q) :: n1 = shape(q,0)
    mdl_ptr = transfer(mdl, mdl_ptr)
    call add_inflow_condition(mdl=mdl_ptr%p, iseg=iseg, coords=coords, t=t, q=q)
end subroutine f90wrap_add_inflow_condition

subroutine f90wrap_set_scheme(mdl, scheme)
    use m_sw_mono, only: model, set_scheme
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    character(*), intent(in), optional :: scheme
    mdl_ptr = transfer(mdl, mdl_ptr)
    call set_scheme(mdl=mdl_ptr%p, scheme=scheme)
end subroutine f90wrap_set_scheme

subroutine f90wrap_model_finalise(mdl)
    use m_sw_mono, only: model, model_finalise
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call model_finalise(mdl=mdl_ptr%p)
    deallocate(mdl_ptr%p)
end subroutine f90wrap_model_finalise

! End of module m_sw_mono defined in file m_sw_mono.f90

