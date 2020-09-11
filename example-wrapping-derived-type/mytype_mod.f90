module mytype_mod

  implicit none

  type mytype
     integer :: a
     character(16) :: strickler_type = "constant"
     real(8), dimension(:), allocatable  ::  m
  end type mytype
  
  type mytype2
     type(mytype), dimension(:), allocatable :: typ
  end type mytype2

contains

  subroutine mytype_initialise(obj,dim_)
    implicit none
    integer, intent(in) :: dim_
    type(mytype), intent(inout) :: obj
    obj%a = dim_
    allocate(obj%m(dim_))
  end subroutine

  subroutine mytype2_initialise(objj,dimens,dim__)
    implicit none
    integer, intent(in) :: dimens
    integer, intent(in) :: dim__
    type(mytype2), intent(out) :: objj
    integer :: i
    allocate(objj%typ(dimens))
    do i = 1, dimens
        call mytype_initialise(objj%typ(i),dim__)
    end do
  end subroutine

  subroutine set_m(obj,m)
    implicit none
    type(mytype), intent(inout) :: obj
    real(8), dimension(:), intent(in) :: m
    if (allocated(obj%m)) deallocate(obj%m)
    obj%a = size(m)
    allocate(obj%m(size(m)))
    obj%m(:) = m
   end subroutine

end module mytype_mod
