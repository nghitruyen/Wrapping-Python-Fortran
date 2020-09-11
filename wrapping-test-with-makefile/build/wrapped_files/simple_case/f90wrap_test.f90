! Module module_test defined in file test.f90

subroutine f90wrap_real_array__array__item(this, nd, dtype, dshape, dloc)
    use module_test, only: real_array
    implicit none
    type real_array_ptr_type
        type(real_array), pointer :: p => NULL()
    end type real_array_ptr_type
    integer, intent(in) :: this(2)
    type(real_array_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%item)
    dloc = loc(this_ptr%p%item)
end subroutine f90wrap_real_array__array__item

subroutine f90wrap_real_array_initialise(this)
    use module_test, only: real_array
    implicit none
    
    type real_array_ptr_type
        type(real_array), pointer :: p => NULL()
    end type real_array_ptr_type
    type(real_array_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_real_array_initialise

subroutine f90wrap_real_array_finalise(this)
    use module_test, only: real_array
    implicit none
    
    type real_array_ptr_type
        type(real_array), pointer :: p => NULL()
    end type real_array_ptr_type
    type(real_array_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_real_array_finalise

subroutine f90wrap_testf(x)
    use module_test, only: testf, real_array
    implicit none
    
    type real_array_ptr_type
        type(real_array), pointer :: p => NULL()
    end type real_array_ptr_type
    type(real_array_ptr_type) :: x_ptr
    integer, intent(in), dimension(2) :: x
    x_ptr = transfer(x, x_ptr)
    call testf(x=x_ptr%p)
end subroutine f90wrap_testf

! End of module module_test defined in file test.f90

