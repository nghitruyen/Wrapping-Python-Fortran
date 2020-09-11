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
!  Fill Ghost Cells with appropriate values depending of the type of the boundary considered
!    ( with eventual pre calculation using Subroutine set_bc)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE calc_boundary_state( mesh , hL , zL , uL , vL , hR , zR , uR , vR )

   USE m_model
 !  USE m_mpi
   USE m_user_data

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   real(rp), intent(in )  ::  hL , uL , vL , zL , zR
   real(rp), intent(out)  ::  hR , uR , vR

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   ib = mesh%edge(ie)%lim

   select case( mesh%edgeb(ib)%typlim )

   !===================================================================================================================!
   !  Neumann Boundary Condition
   !===================================================================================================================!

      case( 'neumann' )

         hR  =  hL

         uR  =  uL

         vR  =  vL

   !===================================================================================================================!
   !  Wall Boundary
   !===================================================================================================================!

      case( 'wall' )

         hR  =   hL + zL - zR

         uR  = - uL

         vR  =   vL

   !===================================================================================================================!
   !  Inflow Boundary  :  'discharg1'
   !===================================================================================================================!

      case( 'discharg1' )

         hR  =  hL

         uR  =  bc%inflow(ib)

         vR  =  vL

   !===================================================================================================================!
   !  Inflow Boundary  :  'discharg2'
   !===================================================================================================================!

      case( 'discharg2' )

         hR  =  bc%inflow(ib)

         uR  =  uL  +  two * sqrt(g) * ( sqrt( hL ) - sqrt( hR ) )

         vR  =  vL

   !===================================================================================================================!
   !  Outflow Boundary  :  'transm'
   !===================================================================================================================!

      case( 'transm' )

         if      ( uL > sqrt( g * hL ) ) then

            hR  =  hL
            uR  =  uL

         else if ( uL > 0._rp .and. hL > heps ) then

            hR  =  hL
            uR  =  uL

         else

            hR  =   hL
            uR  = - uL

         end if

         vR  =  vL

   !===================================================================================================================!
   !  Outflow Boundary  :  'ratcurve'
   !===================================================================================================================!

      case( 'ratcurve' )

         hR  =  bc%outflow(            ib )

         uR  =  bc%outflow( mesh%neb + ib )

         vR  =  vL

   !===================================================================================================================!
   !  Outflow Boundary  :  'zspresc'
   !===================================================================================================================!

      case( 'zspresc' )

         hR  =  max( 0._rp , outflow_user( tc , mesh%cellb(ib)%grav%x , mesh%cellb(ib)%grav%y ) - &
                             bathy_cell( mesh%edge(ie)%cell(2) ) )

         uR  =  uL  +  two * sqrt(g) * ( sqrt( hL ) - sqrt( hR ) )

         vR  =  vL

   !===================================================================================================================!
   !  Outflow Boundary  :  'hpresc'
   !===================================================================================================================!

      case( 'hpresc' )

         hR  =  outflow_user( tc , mesh%cellb(ib)%grav%x , mesh%cellb(ib)%grav%y )

         uR  =  uL  +  two * sqrt(g) * ( sqrt( hL ) - sqrt( hR ) )

         vR  =  vL

   !===================================================================================================================!
   !  Wall Instead
   !===================================================================================================================!

      case default

         hR  =   hL + zL - zR

         uR  = - uL

         vR  =   vL

   end select

END SUBROUTINE calc_boundary_state


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Fill Ghost Unknows in dof structure using above boundary_state Subroutine
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE fill_bc( dof , mesh )

   USE m_model
 !  USE m_mpi
   USE m_user_data

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type( msh ), intent(in   )  ::  mesh
   type( unk ), intent(inout)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp)  ::  hL , uL , vL , zL
   real(rp)  ::  hR , uR , vR , zR

   integer(ip)  ::  iL , iR

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   do ie = 1,mesh%ne

      if ( mesh%edge(ie)%boundary ) then

         iL  =  mesh%edge(ie)%cell(1)
         iR  =  mesh%edge(ie)%cell(2)

         hL  =  dof%h( iL )
         hR  =  dof%h( iR )

         zL  =  bathy_cell( iL )
         zR  =  bathy_cell( iR )

         uL  =  mesh%edge(ie)%normal%x * dof%u( iL ) + &
                mesh%edge(ie)%normal%y * dof%v( iL )

         vL  =  mesh%edge(ie)%normal%x * dof%v( iL ) - &
                mesh%edge(ie)%normal%y * dof%u( iL )

         call calc_boundary_state( mesh , hL , zL , uL , vL , hR , zR , uR , vR )

         dof%h( iR )  =  hR

         dof%u( iR )  =  mesh%edge(ie)%normal%x * uR - mesh%edge(ie)%normal%y * vR

         dof%v( iR )  =  mesh%edge(ie)%normal%y * uR + mesh%edge(ie)%normal%x * vR

      end if

   end do

END SUBROUTINE fill_bc


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Eventual Calculation of Ghost Values
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE set_bc( dof , mesh )

   USE m_numeric
   USE m_model
 !  USE m_mpi
   USE m_user_data

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type( msh ), intent(in)  ::  mesh
   type( unk ), intent(in)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp)  ::  sum_pow_h , zin , qin , qout

   integer(ip)  ::  num_bc

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!


   bc%inflow   =  zero
   bc%outflow  =  zero

   do num_bc = 1,bc%nb

      !================================================================================================================!
      !  discharg1
      !================================================================================================================!

      if ( bc%typ(num_bc,1) == 'discharg1' ) then

         !=============================================================================================================!
         !  loading Qin
         !=============================================================================================================!

         if ( bc%typ(num_bc,2) == 'file' ) then

            qin = linear_interp( bc%hyd( bc%grpf( num_bc ) )%t , &
                                 bc%hyd( bc%grpf( num_bc ) )%q , tc )

         else

            qin = inflow_user( tc , zero , zero )

         end if

         !=============================================================================================================!
         !
         !=============================================================================================================!

         sum_pow_h  =  zero

         do ib = 1,mesh%neb

            if ( mesh%edgeb(ib)%typlim == 'discharg1' .and. mesh%edgeb(ib)%group  == num_bc ) then

               ie  =  mesh%edgeb(ib)%ind

               i  =  mesh%edge(ie)%cell(1)

               if ( dof%h(i) > heps ) sum_pow_h  =  sum_pow_h  +  dof%h(i)**d5p3 * mesh%edge(ie)%length

             end if

         end do

      !   call mpi_sum_r( sum_pow_h )

         !=============================================================================================================!
         !
         !=============================================================================================================!

         do ib = 1,mesh%neb

            if ( mesh%edgeb(ib)%typlim == 'discharg1' .and. mesh%edgeb(ib)%group  == num_bc ) then

               ie  =  mesh%edgeb(ib)%ind

               i  =  mesh%edge(ie)%cell(1)

               if ( dof%h(i) > heps ) then

                  bc%inflow(            ib )  =  -  qin * dof%h(i)**d2p3 / sum_pow_h

                  bc%inflow( mesh%neb + ib )  =     dof%h(i) * bc%inflow( ib )

               end if

            end if

         end do

      end if


      !================================================================================================================!
      !  discharg2
      !================================================================================================================!

      if ( bc%typ(num_bc,1) == 'discharg2' ) then

         !=============================================================================================================!
         !  loading Qin
         !=============================================================================================================!

         if ( bc%typ(num_bc,2) == 'file' ) then

            qin = linear_interp( bc%hyd( bc%grpf( num_bc ) )%t , &
                                 bc%hyd( bc%grpf( num_bc ) )%q , tc )

         else

            qin = inflow_user( tc , zero , zero )

         end if

         !=============================================================================================================!
         !
         !=============================================================================================================!

         call Newton_Qin( qin , dof , mesh , zin )

         !=============================================================================================================!
         !
         !=============================================================================================================!

         do ib = 1,mesh%neb

            if ( mesh%edgeb(ib)%typlim == 'discharg2' .and. mesh%edgeb(ib)%group  == num_bc ) then

               ie  =  mesh%edgeb(ib)%ind

               i  =  mesh%edge(ie)%cell(1)
               j  =  mesh%edge(ie)%cell(2)

               if ( dof%h(i) > heps ) bc%inflow(ib)  =  max( 0._rp , zin  -  bathy_cell(j) )

            end if

         end do

      end if


      !================================================================================================================!
      !  ratcurve
      !================================================================================================================!

      if ( bc%typ(num_bc,1) == 'ratcurve' ) then

         !=============================================================================================================!
         !
         !=============================================================================================================!

         if ( nt == 0 ) then

            qout = zero

            do ib = 1,mesh%neb

               if ( mesh%edgeb(ib)%typlim == 'ratcurve' .and. mesh%edgeb(ib)%group  == num_bc ) then

                  ie  =  mesh%edgeb(ib)%ind

                  i  =  mesh%edge(ie)%cell(1)

                  qout  =  qout  +  dof%h(i) * ( dof%u(i) * mesh%edge(ie)%normal%x + &
                                                 dof%v(i) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length

                end if

            end do

     !       call mpi_sum_r( qout )

         !   bc%rat( bc%grpf( num_bc ) )%zout  =  &

       !     max( 0._rp , linear_interp( bc%rat( bc%grpf( num_bc ) )%q , &
        !                                bc%rat( bc%grpf( num_bc ) )%h , qout ) )

         else

            !TODO relaxation needed, fix coef !?

            bc%rat( bc%grpf( num_bc ) )%zout  =  &

            0.95_rp * bc%rat( bc%grpf( num_bc ) )%zout + &
            0.05_rp * max( 0._rp , linear_interp( bc%rat( bc%grpf( num_bc ) )%q , &
                                                  bc%rat( bc%grpf( num_bc ) )%h , bc%sum_mass_flux( num_bc ) ) )

         end if

         !=============================================================================================================!
         !
         !=============================================================================================================!

         sum_pow_h  =  zero

         do ib = 1,mesh%neb

            if ( mesh%edgeb(ib)%typlim == 'ratcurve' .and. mesh%edgeb(ib)%group  == num_bc ) then

               ie  =  mesh%edgeb(ib)%ind

               i  =  mesh%edge(ie)%cell(1)
               j  =  mesh%edge(ie)%cell(2)

               bc%outflow( ib )  =  max( 0._rp , bc%rat( bc%grpf( num_bc ) )%zout - bathy_cell(j) + &   !TODO
                                                 bc%rat( bc%grpf( num_bc ) )%z_rat_ref )                !right bathy !?

               sum_pow_h  =  sum_pow_h  +  bc%outflow(ib)**d5p3 * mesh%edge(ie)%length

             end if

         end do

  !       call mpi_sum_r( sum_pow_h )

         if ( sum_pow_h > zerom ) then

            do ib = 1,mesh%neb

               if ( mesh%edgeb(ib)%typlim == 'ratcurve' .and. mesh%edgeb(ib)%group  == num_bc ) then

                  ie  =  mesh%edgeb(ib)%ind

                  i  =  mesh%edge(ie)%cell(1)
                  j  =  mesh%edge(ie)%cell(2)

                  qout  =  linear_interp( bc%rat( bc%grpf( num_bc ) )%h(:) , &
                                          bc%rat( bc%grpf( num_bc ) )%q(:) , &
                                          bc%rat( bc%grpf( num_bc ) )%zout )

                  bc%outflow( mesh%neb + ib )  =  bc%outflow(ib)**d2p3 * qout / sum_pow_h

                end if

            end do

         end if

      end if

   end do

END SUBROUTINE set_bc


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE Newton_Qin( qin , dof , mesh , zs )

   USE m_model
 !  USE m_mpi

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type( msh ), intent(in)  ::  mesh
   type( unk ), intent(in)  ::  dof

   real(rp), intent(in )  ::  qin
   real(rp), intent(out)  ::  zs

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   integer(ip)  ::  av

   real(rp)  ::  r , c , z , dzs , s1 , s2

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   c   =  two * sqrt( g )

   zs  =  zero

   av  =  0

   do ie = 1,mesh%neb

      if ( mesh%edgeb(ie)%typlim(1:8) == 'discharg' ) then

         i  =  mesh%edge( mesh%edgeb(ie)%ind )%cell(1)

         if ( dof%h(i) > heps ) then

            av  =  av  +  1

            zs  =  zs  +  bathy_cell(i)  +  dof%h(i)

         end if

      end if

   end do

  ! call mpi_sum_r( zs )
 !  call mpi_sum_i( av )

   if ( av > 0 ) then

      zs  =  zs / real(av,rp)

   else

      zs  =  zero

   end if

   dzs =  one

   k  =  0

   do while ( k <= 50 .and. dzs > 10._rp * zerom )

      dzs  =  zs

      s1  =  zero
      s2  =  zero

      do ie = 1,mesh%neb

         if ( mesh%edgeb(ie)%typlim(1:8) == 'discharg' ) then

            i  =  mesh%edge( mesh%edgeb(ie)%ind )%cell(1)

            if ( dof%h(i) > heps ) then

               r  =  dof%u(i) * mesh%edge( mesh%edgeb(ie)%ind )%normal%x + &
                     dof%v(i) * mesh%edge( mesh%edgeb(ie)%ind )%normal%y + c * sqrt( dof%h(i) )

               z  =  max( zero , zs - bathy_cell(i) )

               s1  =  s1  -  ( z * ( r -        c * sqrt( z ) ) ) * mesh%edge( mesh%edgeb(ie)%ind )%length
               s2  =  s2  -  (       r - d3p2 * c * sqrt( z )   ) * mesh%edge( mesh%edgeb(ie)%ind )%length

            end if

         end if

      end do

   !   call mpi_sum_r( s1 )
    !  call mpi_sum_r( s2 )

      zs = zs - ( s1 - qin ) / s2

      dzs = abs( 1._rp - dzs / zs )

      k = k + 1

   end do

 !  if ( dzs > 10._rp * zerom ) call Stopping_Program_Sub( 'Problem of convergence in Newton_Qin' )

END SUBROUTINE Newton_Qin
