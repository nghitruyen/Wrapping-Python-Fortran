subroutine f90wrap_foo(a, b)
    implicit none
    external foo
    
    real(4), intent(in) :: a
    integer :: b
    call foo(a, b)
end subroutine f90wrap_foo

