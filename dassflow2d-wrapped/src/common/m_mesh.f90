!======================================================================================================================!
!
!                    DassFlow Version 2.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures:
!    - Shallow Module (Shallow Water Model, Finite Volume Method), i.e. the present code.
!    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current main developer:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!  with scientific and/or programming contributions of:
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!               K. Larnier (Fluid Mechanics Institute of Toulouse IMFT).
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!  and former other developers (M. Honnorat and J. Marin).
!
!  Scientific Contact : jerome.monnier@insa-toulouse.fr
!  Technical  Contact : frederic.couderc@math.univ-toulouse.fr
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_mesh
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_mesh

   USE m_common
   USE m_linear_algebra

   implicit none

   !===================================================================================================================!
   !  Maximum number of Nodes/Edges connecting the Cell
   !===================================================================================================================!

   integer, parameter  ::  maxed = 4               ! Maximum number of Nodes/Edges connecting the Cell

   !===================================================================================================================!
   !  Node Structure
   !===================================================================================================================!

   TYPE NodeType

      integer(ip), dimension(:), allocatable  ::  cell      ! Connected Cells Indexes

      integer(ip), dimension(:), allocatable  ::  edge      ! Connected Edges Indexes

      type(vec2d)  ::  coord                                ! Node coordinates

      logical  ::  boundary                                 ! true boolean if Node is a at Mesh Boundary

      integer(ip)  ::  lim                                  ! Boundary Node Index

   END TYPE NodeType

   !===================================================================================================================!
   !  Boundary Node Structure
   !===================================================================================================================!

   TYPE NodeTypeLim

      integer(ip)  ::  ind                         ! Global Node Index

      character(len=lchar)  ::  typlim             ! Type of the Boundary condition

      integer(ip)  ::  group                       ! Group Number in case of Multi Boundary Conditions

   END TYPE NodeTypeLim

   !===================================================================================================================!
   !  Cell Structure
   !===================================================================================================================!

   TYPE CellType

      integer(ip), dimension(maxed)  ::  node      ! Suited list of Nodes connecting the Cell

      integer(ip), dimension(maxed)  ::  cell      ! Neighboring Cells Index (referenced to upper Nodes list)

      integer(ip), dimension(maxed)  ::  edge      ! Edges Indexes

      integer(ip)  ::  nbed                        ! Number of Neighboring Cells/Edges

      logical  ::  boundary                        ! True boolean if Cell is a at Mesh Boundary

      real(rp)  ::  surf                           ! Cell surface

      real(rp)  ::  invsurf                        ! Inverse of Cell surface

      real(rp)  ::  peri                           ! Cell perimeter

      type(vec2d)  ::  grav                        ! Cell gravity center

   END TYPE CellType

   !===================================================================================================================!
   !  Boundary Cell Structure
   !===================================================================================================================!

   TYPE CellTypeLim

      integer(ip)  ::  ind                         ! Global Cell Index

      character(len=lchar)  ::  typlim             ! Type of the Boundary Condition

      integer(ip)  ::  group                       ! Group Number in case of Multi Boundary Conditions

      integer(ip)  ::  cell                        ! Neighboring Cell

      type(vec2d)  ::  grav                        ! Cell gravity center

   END TYPE CellTypeLim

   !===================================================================================================================!
   !  Edge Structure
   !===================================================================================================================!

   TYPE EdgeType

      integer(ip)  ::  node(2)                     ! Node Indexes connecting the Edge

      integer(ip)  ::  cell(2)                     ! Cell Indexes linked to the Edge

      logical  ::  boundary                        ! True boolean if Edge is a at mesh Boundary

      logical  ::  subdomain                       ! True boolean if Edge is a at sub-domain Boundary

      integer(ip)  ::  lim                         ! Boundary Edge Index

      real(rp)  ::  length                         ! Edge length

      type(vec2d)  ::  center                      ! Edge center

      type(vec2d)  ::  normal                      ! Edge normal (oriented from 1 to 2 linked cells)

      type(vec2d)  ::  tangent                     ! Edge tangent (oriented from 1 to 2 linked nodes)

      type(vec2d)  ::  vcell                       ! Cells(1&2) gravity center vector (oriented from 1 to 2)

      type(vec2d)  ::  v_edge_cell(2)              ! Vector going from the edge center to the cell(1&2) gravity center

   END TYPE EdgeType

   !===================================================================================================================!
   !  Boundary Edge Structure
   !===================================================================================================================!

   TYPE EdgeTypeLim

      integer(ip)  ::  ind                         ! Global Edge Index

      character(len=lchar)  ::  typlim             ! Type of the Boundary condition

      integer(ip)  ::  group                       ! Group Number in case of Multi Boundary Conditions

      integer(ip)  ::  perio                       ! Symetric Edge Index in case of Periodic Boundary Condition

   END TYPE EdgeTypeLim

   !===================================================================================================================!
   !  Global Mesh Structure
   !===================================================================================================================!

   TYPE msh

      integer(ip)  ::  nn                                            ! Number of mesh Nodes
      integer(ip)  ::  nnb                                           ! Number of mesh Nodes at Boundary

      integer(ip)  ::  nc                                            ! Number of mesh Cells
      integer(ip)  ::  ncb                                           ! Number of mesh Cells at Boundary

      integer(ip)  ::  ne                                            ! Number of mesh Edges
      integer(ip)  ::  neb                                           ! Number of mesh Edges at Boundary

      type( NodeType    ), dimension(:), allocatable  ::  node       ! Node Structure
      type( NodeTypeLim ), dimension(:), allocatable  ::  nodeb      ! Boundary Node Structure

      type( CellType    ), dimension(:), allocatable  ::  cell       ! Cell Structure
      type( CellTypeLim ), dimension(:), allocatable  ::  cellb      ! Boundary Cell Structure

      type( EdgeType    ), dimension(:), allocatable  ::  edge       ! Edge Structure
      type( EdgeTypeLim ), dimension(:), allocatable  ::  edgeb      ! Boundary Edge Structure

      character(len=lchar)  ::  file_name                            ! Name of the mesh file

      real(rp)  ::  scal                                             ! Mesh scaling factor

      real(rp)  ::  surf                                             ! Mesh total surface

   END TYPE msh

   !===================================================================================================================!
   !  Point in Mesh Structure
   !===================================================================================================================!

   TYPE point_in_mesh

      type( vec2d )  ::  coord

      integer(ip)  ::  cell

   END TYPE point_in_mesh

CONTAINS

!===================================================
!	
!===================================================
   subroutine NodeType_initialise(node,dim1,dim2)
      implicit none
      type(NodeType), intent(out) :: node
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
      allocate(node%cell(dim1))
      allocate(node%edge(dim2))
      node%boundary = .true.
      node%lim = 0
      node%coord%x = 0
      node%coord%y = 0
   end subroutine

   subroutine NodeTypeLim_initialise(nodelim,ind,group)
      implicit none
      integer, intent(in) :: ind
      integer, intent(in) :: group
      type(NodeTypeLim), intent(out) :: nodelim
      nodelim%ind = ind
      nodelim%group = group
      nodelim%typlim = 'NULL'
   end subroutine

   subroutine msh_initialise(mesh)
      implicit none
      type(msh), intent(out) :: mesh
      mesh%nn   =  nx * ny
      mesh%ne   =  nx * (ny-1) + ny * (nx-1)
      mesh%nc   =  (nx-1) * (ny-1)
      mesh%nnb  =  2 * ( (nx-2) + (ny-2) ) + 4
      mesh%neb  =  2 * ( (nx-1) + (ny-1) )
      mesh%ncb  =  mesh%neb
      allocate( mesh%node( mesh%nn ) )
      allocate( mesh%edge( mesh%ne ) )
      allocate( mesh%cell( mesh%nc ) )
      allocate( mesh%nodeb( mesh%nnb ) )
      allocate( mesh%edgeb( mesh%neb ) )
      allocate( mesh%cellb( mesh%ncb ) )
      dx  =  lx / float( nx - 1 )
      dy  =  ly / float( ny - 1 )
      k  = 0
      kb = 0
      do j = 1,ny
         do i = 1,nx
            k = k + 1
            mesh%node(k)%coord%x  =  float( i - 1 ) * dx
            mesh%node(k)%coord%y  =  float( j - 1 ) * dy
            if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
               kb = kb + 1
               mesh%node(k)%boundary  =  .true.
               mesh%node(k)%lim  =  kb
               mesh%nodeb(kb)%ind  =  k
            else
               mesh%node(k)%boundary  =  .false.
               mesh%node(k)%lim  =  0
            end if
         end do
      end do

      k  = 0
      kb = 0

      do j = 1,ny-1
         do i = 1,nx

            k = k + 1

            mesh%edge(k)%node(1)  =  i + (j-1) * nx
            mesh%edge(k)%node(2)  =  i +  j    * nx

            if      ( i == 1 ) then

               kb = kb + 1
  
               mesh%edge(k)%boundary  =  .true.

               mesh%edge(k)%lim  =  kb

               mesh%edgeb(kb)%ind  =  k

               mesh%edgeb(kb)%typlim  =  bc_W

               mesh%edgeb(kb)%group  =  1

               mesh%edge(k)%cell(1)  =  i + (j-1) * (nx-1)

               mesh%edge(k)%cell(2)  =  mesh%nc + kb

               if ( bc_W == 'periodic' ) then

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) + nx - 2

                  mesh%edgeb(kb)%perio = k + nx - 1

               else

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)

                  mesh%edgeb(kb)%perio = -1

               end if

            else if ( i == nx ) then

               kb = kb + 1

               mesh%edge(k)%boundary  =  .true.

               mesh%edge(k)%lim  =  kb

               mesh%edgeb(kb)%ind  =  k

               mesh%edgeb(kb)%typlim  =  bc_E

               mesh%edgeb(kb)%group  =  2

               mesh%edge(k)%cell(1)  =  i-1 + (j-1) * (nx-1)

               mesh%edge(k)%cell(2)  =  mesh%nc + kb

               if ( bc_E == 'periodic' ) then

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) - nx + 2

                  mesh%edgeb(kb)%perio = k - nx + 1
 
               else

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)

                  mesh%edgeb(kb)%perio = -1

               end if

            else

               mesh%edge(k)%boundary  =  .false.

               mesh%edge(k)%lim  =  0

               mesh%edge(k)%cell(1)  =  i-1 + (j-1) * (nx-1)
               mesh%edge(k)%cell(2)  =  i   + (j-1) * (nx-1)

            end if

         end do
      end do

      do i = 1,nx-1
         do j = 1,ny

            k = k + 1

            mesh%edge(k)%node(1)  =  i   + (j-1) * nx
            mesh%edge(k)%node(2)  =  i+1 + (j-1) * nx

            if      ( j == 1 ) then

               kb = kb + 1

               mesh%edge(k)%boundary  =  .true.

               mesh%edge(k)%lim  =  kb

               mesh%edgeb(kb)%ind  =  k

               mesh%edgeb(kb)%typlim  =  bc_S

               mesh%edgeb(kb)%group  =  3

               mesh%edge(k)%cell(1)  =  i + j - 1

               mesh%edge(k)%cell(2)  =  mesh%nc + kb

               if ( bc_S == 'periodic' ) then

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) + ( nx - 1 ) * ( ny - 2 )

                  mesh%edgeb(kb)%perio = k + ny - 1

               else

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)

                  mesh%edgeb(kb)%perio = -1

               end if

            else if ( j == ny ) then

               kb = kb + 1

               mesh%edge(k)%boundary  =  .true.

               mesh%edge(k)%lim  =  kb

               mesh%edgeb(kb)%ind  =  k

               mesh%edgeb(kb)%typlim  =  bc_N

               mesh%edgeb(kb)%group  =  4

               mesh%edge(k)%cell(1)  =  i + (j-2) * (nx-1)

               mesh%edge(k)%cell(2)  =  mesh%nc + kb

               if ( bc_N == 'periodic' ) then

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) - ( nx - 1 ) * ( ny - 2 )

                  mesh%edgeb(kb)%perio = k - ny + 1

               else

                  mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)

                  mesh%edgeb(kb)%perio = -1

               end if

            else

               mesh%edge(k)%boundary  =  .false.

               mesh%edge(k)%lim  =  0

               mesh%edge(k)%cell(1)  =  i + (j-2) * (nx-1)
               mesh%edge(k)%cell(2)  =  i + (j-1) * (nx-1)

            end if

         end do
      end do

   !===================================================================================================================!
   !  Cells variables (boundary excluded from global indexes, see below)
   !===================================================================================================================!

      k  = 0

      do j = 1,ny-1
         do i = 1,nx-1

            k = k + 1

            mesh%cell(k)%node(1)  =  k + j - 1
            mesh%cell(k)%node(2)  =  k + j
            mesh%cell(k)%node(3)  =  k + nx + j
            mesh%cell(k)%node(4)  =  k + nx + j - 1

            mesh%cell(k)%cell(1)  =  k - (nx-1)
            mesh%cell(k)%cell(2)  =  k + 1
            mesh%cell(k)%cell(3)  =  k + (nx-1)
            mesh%cell(k)%cell(4)  =  k - 1

            mesh%cell(k)%edge(1)  =  j + nx * ( ny - 1 ) + ny * ( i - 1 )
            mesh%cell(k)%edge(2)  =  k + j
            mesh%cell(k)%edge(3)  =  mesh%cell(k)%edge(1) + 1
            mesh%cell(k)%edge(4)  =  mesh%cell(k)%edge(2) - 1

            if ( i == 1 .or. i == nx-1 .or. j == 1 .or. j == ny-1 ) then

               mesh%cell(k)%boundary  =  .true.

               if ( j == 1    ) mesh%cell(k)%cell(1)  =  mesh%edge( mesh%cell(k)%edge(1) )%cell(2)
               if ( i == nx-1 ) mesh%cell(k)%cell(2)  =  mesh%edge( mesh%cell(k)%edge(2) )%cell(2)
               if ( j == ny-1 ) mesh%cell(k)%cell(3)  =  mesh%edge( mesh%cell(k)%edge(3) )%cell(2)
               if ( i == 1    ) mesh%cell(k)%cell(4)  =  mesh%edge( mesh%cell(k)%edge(4) )%cell(2)

            else

               mesh%cell(k)%boundary  =  .false.

            end if

            mesh%cell(k)%nbed  =  4

         end do
      end do
   end subroutine

!=============================================================================
!		TO VISUALIZE FOR IMPORTING DATA
!==============================================================================

   subroutine display_mesh_cell(mesh)
      implicit none
      type(msh), intent(in) :: mesh
      write(*,*) mesh%cell(i)%invsurf
      write(*,*) mesh%cell(i)%surf
      write(*,*) mesh%cell(i)%boundary
      do i=1,mesh%nc
         !write(*,*) mesh%cell(i)%grav%x
        ! write(*,*) mesh%cell(i)%grav%y
      end do
    end subroutine
   !===================================================================================================================!
   !  Nodes variables (boundary included in global indexes)
   !===================================================================================================================!

  
!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of one Cell variable
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE alloc_cell( var , mesh )

      implicit none

      type( msh ), intent(in)  ::  mesh

      real(rp), dimension(:), allocatable, intent(out)  ::  var

      allocate( var( mesh%nc ) )

      var(:)  =  0._rp

   END SUBROUTINE alloc_cell


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of one Edge variable
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE alloc_edge( var , mesh )

      implicit none

      type( msh ), intent(in)  ::  mesh

      real(rp), dimension(:), allocatable, intent(out)  ::  var

      allocate( var( mesh%ne ) )

      var(:)  =  0._rp

   END SUBROUTINE alloc_edge


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of one Node variable
!
!!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE alloc_node( var , mesh )

      implicit none

      type( msh ), intent(in)  ::  mesh

      real(rp), dimension(:), allocatable, intent(out)  ::  var

      allocate( var( mesh%nn ) )

      var(:)  =  0._rp

   END SUBROUTINE alloc_node


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Increase or Decrease Allocated Mesh Array Memory
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE reallocate_cell( var , new )

      implicit none

      type( CellType ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( CellType ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_cell


   SUBROUTINE reallocate_edge( var , new )

      implicit none

      type( EdgeType ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( EdgeType ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_edge


   SUBROUTINE reallocate_node( var , new )

      implicit none

      type( NodeType ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( NodeType ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_node


   SUBROUTINE reallocate_cellb( var , new )

      implicit none

      type( CellTypeLim ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( CellTypeLim ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_cellb


   SUBROUTINE reallocate_edgeb( var , new )

      implicit none

      type( EdgeTypeLim ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( EdgeTypeLim ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_edgeb


   SUBROUTINE reallocate_nodeb( var , new )

      implicit none

      type( NodeTypeLim ), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      type( NodeTypeLim ), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( : )

         call move_alloc( temp , var )

      end if

   END SUBROUTINE reallocate_nodeb


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of mesh structure
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE dealloc_mesh( mesh )

      implicit none

      type( msh ), intent(inout)  ::  mesh

      if ( allocated( mesh%cell  ) ) deallocate( mesh%cell  )
      if ( allocated( mesh%cellb ) ) deallocate( mesh%cellb )
      if ( allocated( mesh%edge  ) ) deallocate( mesh%edge  )
      if ( allocated( mesh%edgeb ) ) deallocate( mesh%edgeb )
      if ( allocated( mesh%node  ) ) deallocate( mesh%node  )
      if ( allocated( mesh%nodeb ) ) deallocate( mesh%nodeb )

   END SUBROUTINE dealloc_mesh


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calculation of cells connectivity
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE calc_cells_connectivity( mesh )

      implicit none

      type( msh ), intent(inout)  ::  mesh

      integer(ip)  ::  side1 , side2 , tri1 , tri2 , col( 4 , maxed * mesh%nc ) , node( maxed ) , n1 , n2

      mesh%cell(:)%nbed = maxed

      j = 0

      do i = 1,mesh%nc

         mesh%cell(i)%cell(:) = 0

         node(:) = mesh%cell(i)%node(:)

         do k = 1,maxed

            n1 =     k
            n2 = mod(k,maxed) + 1

            if ( node(n1) == node(n2) ) then

               mesh%cell(i)%cell(k) = -1

               mesh%cell(i)%nbed = mesh%cell(i)%nbed - 1

            else if ( node(n1) < node(n2) ) then

               j = j + 1 ; col(1:4,j) = (/ node(n1), node(n2), k , i /)

            else

               j = j + 1 ; col(1:4,j) = (/ node(n2), node(n1), k , i /)

            end if

         end do

      end do

      call i4col_sort_a( 4 , j , col )

      mesh%ne  = 0
      mesh%neb = 0

      i = 1

      do

         if ( j <= i ) then

            if ( i == j ) mesh%neb = mesh%neb + 1

            exit

         end if

         mesh%ne = mesh%ne + 1

         if ( col(1,i) /= col(1,i+1) .or. col(2,i) /= col(2,i+1) ) then

            i = i + 1

            mesh%neb = mesh%neb + 1

            cycle

         end if

         side1  =  col(3,i  )
         side2  =  col(3,i+1)

         tri1   =  col(4,i  )
         tri2   =  col(4,i+1)

         mesh%cell( tri1 )%cell( side1 )  =  tri2
         mesh%cell( tri2 )%cell( side2 )  =  tri1

         i = i + 2

      end do

   END SUBROUTINE calc_cells_connectivity


END MODULE m_mesh
