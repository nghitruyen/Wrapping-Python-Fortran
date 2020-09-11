! Module dual_num_auto_diff defined in file DNAD.f90

subroutine f90wrap_dual_num__get__x_ad_(this, f90wrap_x_ad_)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in)   :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_x_ad_
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x_ad_ = this_ptr%p%x_ad_
end subroutine f90wrap_dual_num__get__x_ad_

subroutine f90wrap_dual_num__set__x_ad_(this, f90wrap_x_ad_)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in)   :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_x_ad_
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x_ad_ = f90wrap_x_ad_
end subroutine f90wrap_dual_num__set__x_ad_

subroutine f90wrap_dual_num__array__xp_ad_(this, nd, dtype, dshape, dloc)
    use dual_num_auto_diff, only: dual_num
    implicit none
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    integer, intent(in) :: this(2)
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%xp_ad_)
    dloc = loc(this_ptr%p%xp_ad_)
end subroutine f90wrap_dual_num__array__xp_ad_

subroutine f90wrap_dual_num_initialise(this)
    use dual_num_auto_diff, only: dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_dual_num_initialise

subroutine f90wrap_dual_num_finalise(this)
    use dual_num_auto_diff, only: dual_num
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_dual_num_finalise

subroutine f90wrap_dual_num_auto_diff__get__NDV_AD(f90wrap_NDV_AD)
    use dual_num_auto_diff, only: dual_num_auto_diff_NDV_AD => NDV_AD
    implicit none
    integer(4), intent(out) :: f90wrap_NDV_AD
    
    f90wrap_NDV_AD = dual_num_auto_diff_NDV_AD
end subroutine f90wrap_dual_num_auto_diff__get__NDV_AD

! End of module dual_num_auto_diff defined in file DNAD.f90

