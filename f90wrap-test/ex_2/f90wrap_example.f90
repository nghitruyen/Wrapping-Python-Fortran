! Module class_example defined in file ./example.f90

subroutine f90wrap_example__get__first(this, f90wrap_first)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_first
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_first = this_ptr%p%first
end subroutine f90wrap_example__get__first

subroutine f90wrap_example__set__first(this, f90wrap_first)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_first
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%first = f90wrap_first
end subroutine f90wrap_example__set__first

subroutine f90wrap_example__get__second(this, f90wrap_second)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_second
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_second = this_ptr%p%second
end subroutine f90wrap_example__get__second

subroutine f90wrap_example__set__second(this, f90wrap_second)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_second
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%second = f90wrap_second
end subroutine f90wrap_example__set__second

subroutine f90wrap_example__get__third(this, f90wrap_third)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_third
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_third = this_ptr%p%third
end subroutine f90wrap_example__get__third

subroutine f90wrap_example__set__third(this, f90wrap_third)
    use class_example, only: example
    implicit none
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer, intent(in)   :: this(2)
    type(example_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_third
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%third = f90wrap_third
end subroutine f90wrap_example__set__third

subroutine f90wrap_example_initialise(this)
    use class_example, only: example
    implicit none
    
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    type(example_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_example_initialise

subroutine f90wrap_example_finalise(this)
    use class_example, only: example
    implicit none
    
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    type(example_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_example_finalise

subroutine f90wrap_return_example_first(ret_instance, first)
    use class_example, only: return_example, example
    implicit none
    
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    type(example_ptr_type) :: ret_instance_ptr
    integer, intent(out), dimension(2) :: ret_instance
    integer :: first
    allocate(ret_instance_ptr%p)
    ret_instance_ptr%p = return_example(first=first)
    ret_instance = transfer(ret_instance_ptr, ret_instance)
end subroutine f90wrap_return_example_first

subroutine f90wrap_return_example_second(first, ret_instance, second)
    use class_example, only: return_example, example
    implicit none
    
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer :: first
    type(example_ptr_type) :: ret_instance_ptr
    integer, intent(out), dimension(2) :: ret_instance
    integer :: second
    allocate(ret_instance_ptr%p)
    ret_instance_ptr%p = return_example(first=first, second=second)
    ret_instance = transfer(ret_instance_ptr, ret_instance)
end subroutine f90wrap_return_example_second

subroutine f90wrap_return_example_third(first, second, ret_instance, third)
    use class_example, only: return_example, example
    implicit none
    
    type example_ptr_type
        type(example), pointer :: p => NULL()
    end type example_ptr_type
    integer :: first
    integer :: second
    type(example_ptr_type) :: ret_instance_ptr
    integer, intent(out), dimension(2) :: ret_instance
    integer :: third
    allocate(ret_instance_ptr%p)
    ret_instance_ptr%p = return_example(first=first, second=second, third=third)
    ret_instance = transfer(ret_instance_ptr, ret_instance)
end subroutine f90wrap_return_example_third

! End of module class_example defined in file ./example.f90

