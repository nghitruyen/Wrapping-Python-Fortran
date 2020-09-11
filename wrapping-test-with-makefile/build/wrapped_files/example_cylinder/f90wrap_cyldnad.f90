! Module mcyldnad defined in file cyldnad.f90

subroutine f90wrap_cyldnad(vol, radius, height)
    use dual_num_auto_diff, only: dual_num
    use mcyldnad, only: cyldnad
    implicit none
    
    type dual_num_ptr_type
        type(dual_num), pointer :: p => NULL()
    end type dual_num_ptr_type
    type(dual_num_ptr_type) :: vol_ptr
    integer, intent(out), dimension(2) :: vol
    type(dual_num_ptr_type) :: radius_ptr
    integer, intent(in), dimension(2) :: radius
    type(dual_num_ptr_type) :: height_ptr
    integer, intent(in), dimension(2) :: height
    radius_ptr = transfer(radius, radius_ptr)
    height_ptr = transfer(height, height_ptr)
    allocate(vol_ptr%p)
    call cyldnad(vol=vol_ptr%p, radius=radius_ptr%p, height=height_ptr%p)
    vol = transfer(vol_ptr, vol)
end subroutine f90wrap_cyldnad

! End of module mcyldnad defined in file cyldnad.f90

