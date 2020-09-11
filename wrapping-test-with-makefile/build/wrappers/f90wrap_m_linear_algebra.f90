! Module m_linear_algebra defined in file m_linear_algebra.f90

subroutine f90wrap_vec2d__get__x(this, f90wrap_x)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_vec2d__get__x

subroutine f90wrap_vec2d__set__x(this, f90wrap_x)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x = f90wrap_x
end subroutine f90wrap_vec2d__set__x

subroutine f90wrap_vec2d__get__y(this, f90wrap_y)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_vec2d__get__y

subroutine f90wrap_vec2d__set__y(this, f90wrap_y)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_vec2d__set__y

subroutine f90wrap_vec2d_initialise(this)
    use m_linear_algebra, only: vec2d
    implicit none
    
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    type(vec2d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_vec2d_initialise

subroutine f90wrap_vec2d_finalise(this)
    use m_linear_algebra, only: vec2d
    implicit none
    
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    type(vec2d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_vec2d_finalise

subroutine f90wrap_matrix__array__m(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: matrix
    implicit none
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    integer, intent(in) :: this(2)
    type(matrix_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%m)) then
        dshape(1:2) = shape(this_ptr%p%m)
        dloc = loc(this_ptr%p%m)
    else
        dloc = 0
    end if
end subroutine f90wrap_matrix__array__m

subroutine f90wrap_matrix_initialise(this)
    use m_linear_algebra, only: matrix
    implicit none
    
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    type(matrix_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_matrix_initialise

subroutine f90wrap_matrix_finalise(this)
    use m_linear_algebra, only: matrix
    implicit none
    
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    type(matrix_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_matrix_finalise

subroutine f90wrap_matrixcsr__get__n(this, f90wrap_n)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_matrixcsr__get__n

subroutine f90wrap_matrixcsr__set__n(this, f90wrap_n)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_matrixcsr__set__n

subroutine f90wrap_matrixcsr__get__nnz(this, f90wrap_nnz)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nnz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nnz = this_ptr%p%nnz
end subroutine f90wrap_matrixcsr__get__nnz

subroutine f90wrap_matrixcsr__set__nnz(this, f90wrap_nnz)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nnz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nnz = f90wrap_nnz
end subroutine f90wrap_matrixcsr__set__nnz

subroutine f90wrap_matrixcsr__array__irow(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in) :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%irow)) then
        dshape(1:1) = shape(this_ptr%p%irow)
        dloc = loc(this_ptr%p%irow)
    else
        dloc = 0
    end if
end subroutine f90wrap_matrixcsr__array__irow

subroutine f90wrap_matrixcsr__array__icol(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in) :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%icol)) then
        dshape(1:1) = shape(this_ptr%p%icol)
        dloc = loc(this_ptr%p%icol)
    else
        dloc = 0
    end if
end subroutine f90wrap_matrixcsr__array__icol

subroutine f90wrap_matrixcsr__array__anz(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: matrixcsr
    implicit none
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    integer, intent(in) :: this(2)
    type(matrixcsr_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%anz)) then
        dshape(1:1) = shape(this_ptr%p%anz)
        dloc = loc(this_ptr%p%anz)
    else
        dloc = 0
    end if
end subroutine f90wrap_matrixcsr__array__anz

subroutine f90wrap_matrixcsr_initialise(mat, n, nnz)
    use m_linear_algebra, only: matrixcsr_initialise, matrixcsr
    implicit none
    
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    type(matrixcsr_ptr_type) :: mat_ptr
    integer, intent(out), dimension(2) :: mat
    integer(4), intent(in) :: n
    integer(4), intent(in) :: nnz
    allocate(mat_ptr%p)
    call matrixcsr_initialise(mat=mat_ptr%p, n=n, nnz=nnz)
    mat = transfer(mat_ptr, mat)
end subroutine f90wrap_matrixcsr_initialise

subroutine f90wrap_csr_x_vec_add(a, x, y, n0, n1)
    use m_linear_algebra, only: csr_x_vec_add, matrixcsr
    implicit none
    
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    type(matrixcsr_ptr_type) :: a_ptr
    integer, intent(in), dimension(2) :: a
    real(8), intent(in), dimension(n0) :: x
    real(8), intent(inout), dimension(n1) :: y
    integer :: n0
    !f2py intent(hide), depend(x) :: n0 = shape(x,0)
    integer :: n1
    !f2py intent(hide), depend(y) :: n1 = shape(y,0)
    a_ptr = transfer(a, a_ptr)
    call csr_x_vec_add(A=a_ptr%p, X=x, Y=y)
end subroutine f90wrap_csr_x_vec_add

subroutine f90wrap_matrixcsr_finalise(this)
    use m_linear_algebra, only: matrixcsr
    implicit none
    
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    type(matrixcsr_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_matrixcsr_finalise

subroutine f90wrap_matrixblock__get__n(this, f90wrap_n)
    use m_linear_algebra, only: matrixblock
    implicit none
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixblock_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_matrixblock__get__n

subroutine f90wrap_matrixblock__set__n(this, f90wrap_n)
    use m_linear_algebra, only: matrixblock
    implicit none
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    integer, intent(in)   :: this(2)
    type(matrixblock_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_matrixblock__set__n

subroutine f90wrap_matrixblock__array_getitem__blocks(f90wrap_this, f90wrap_i, &
    blocksitem)
    
    use m_linear_algebra, only: matrix, matrixblock
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(matrixblock_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: blocksitem(2)
    type(matrix_ptr_type) :: blocks_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%blocks)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%blocks)) then
            call f90wrap_abort("array index out of range")
        else
            blocks_ptr%p => this_ptr%p%blocks(f90wrap_i)
            blocksitem = transfer(blocks_ptr,blocksitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_matrixblock__array_getitem__blocks

subroutine f90wrap_matrixblock__array_setitem__blocks(f90wrap_this, f90wrap_i, &
    blocksitem)
    
    use m_linear_algebra, only: matrix, matrixblock
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(matrixblock_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: blocksitem(2)
    type(matrix_ptr_type) :: blocks_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%blocks)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%blocks)) then
            call f90wrap_abort("array index out of range")
        else
            blocks_ptr = transfer(blocksitem,blocks_ptr)
            this_ptr%p%blocks(f90wrap_i) = blocks_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_matrixblock__array_setitem__blocks

subroutine f90wrap_matrixblock__array_len__blocks(f90wrap_this, f90wrap_n)
    
    use m_linear_algebra, only: matrix, matrixblock
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(matrixblock_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%blocks)) then
        f90wrap_n = size(this_ptr%p%blocks)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_matrixblock__array_len__blocks

subroutine f90wrap_matrixblock_initialise(mat, n, sizes, n0)
    use m_linear_algebra, only: matrixblock, matrixblock_initialise
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type(matrixblock_ptr_type) :: mat_ptr
    integer, intent(out), dimension(2) :: mat
    integer(4), intent(in) :: n
    integer(4), intent(in), dimension(n0) :: sizes
    integer :: n0
    !f2py intent(hide), depend(sizes) :: n0 = shape(sizes,0)
    allocate(mat_ptr%p)
    call matrixblock_initialise(mat=mat_ptr%p, n=n, sizes=sizes)
    mat = transfer(mat_ptr, mat)
end subroutine f90wrap_matrixblock_initialise

subroutine f90wrap_matrixblock_copy(src, dst)
    use m_linear_algebra, only: matrixblock_copy, matrixblock
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type(matrixblock_ptr_type) :: src_ptr
    integer, intent(in), dimension(2) :: src
    type(matrixblock_ptr_type) :: dst_ptr
    integer, intent(out), dimension(2) :: dst
    src_ptr = transfer(src, src_ptr)
    allocate(dst_ptr%p)
    call matrixblock_copy(src=src_ptr%p, dst=dst_ptr%p)
    dst = transfer(dst_ptr, dst)
end subroutine f90wrap_matrixblock_copy

subroutine f90wrap_matrixcsr_from_matrixblock(matblock, matcsr, threshold)
    use m_linear_algebra, only: matrixcsr_from_matrixblock, matrixblock, matrixcsr
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    type(matrixblock_ptr_type) :: matblock_ptr
    integer, intent(in), dimension(2) :: matblock
    type(matrixcsr_ptr_type) :: matcsr_ptr
    integer, intent(in), dimension(2) :: matcsr
    real(8), intent(in), optional :: threshold
    matblock_ptr = transfer(matblock, matblock_ptr)
    matcsr_ptr = transfer(matcsr, matcsr_ptr)
    call matrixcsr_from_matrixblock(matblock=matblock_ptr%p, matcsr=matcsr_ptr%p, &
        threshold=threshold)
end subroutine f90wrap_matrixcsr_from_matrixblock

subroutine f90wrap_matrixblock_finalise(this)
    use m_linear_algebra, only: matrixblock
    implicit none
    
    type matrixblock_ptr_type
        type(matrixblock), pointer :: p => NULL()
    end type matrixblock_ptr_type
    type(matrixblock_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_matrixblock_finalise

subroutine f90wrap_matrixcsr_from_numpy_array(array, mat, threshold, n0, n1)
    use m_linear_algebra, only: matrixcsr, matrixcsr_from_numpy_array
    implicit none
    
    type matrixcsr_ptr_type
        type(matrixcsr), pointer :: p => NULL()
    end type matrixcsr_ptr_type
    real(8), intent(in), dimension(n0,n1) :: array
    type(matrixcsr_ptr_type) :: mat_ptr
    integer, intent(out), dimension(2) :: mat
    real(8), intent(in), optional :: threshold
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    allocate(mat_ptr%p)
    call matrixcsr_from_numpy_array(array=array, mat=mat_ptr%p, threshold=threshold)
    mat = transfer(mat_ptr, mat)
end subroutine f90wrap_matrixcsr_from_numpy_array

subroutine f90wrap_cholesky_inplace(n, a, n0, n1)
    use m_linear_algebra, only: cholesky_inplace
    implicit none
    
    integer(4), intent(in) :: n
    real(8), dimension(n0,n1) :: a
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    call cholesky_inplace(n=n, A=a)
end subroutine f90wrap_cholesky_inplace

subroutine f90wrap_solve_using_cholesky_inplace(n, mat, b, x, n0, n1)
    use m_linear_algebra, only: matrix, solve_using_cholesky_inplace
    implicit none
    
    type matrix_ptr_type
        type(matrix), pointer :: p => NULL()
    end type matrix_ptr_type
    integer(4), intent(in) :: n
    type(matrix_ptr_type) :: mat_ptr
    integer, intent(in), dimension(2) :: mat
    real(8), intent(in), dimension(n0) :: b
    real(8), intent(inout), dimension(n1) :: x
    integer :: n0
    !f2py intent(hide), depend(b) :: n0 = shape(b,0)
    integer :: n1
    !f2py intent(hide), depend(x) :: n1 = shape(x,0)
    mat_ptr = transfer(mat, mat_ptr)
    call solve_using_cholesky_inplace(n=n, mat=mat_ptr%p, b=b, x=x)
end subroutine f90wrap_solve_using_cholesky_inplace

! End of module m_linear_algebra defined in file m_linear_algebra.f90

