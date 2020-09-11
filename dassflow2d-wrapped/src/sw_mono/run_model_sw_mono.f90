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
!  Main Routine to run a Shallow-Water simulation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE run_model( mesh , dof0 , dof , cost )

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   !USE m_mpi
  ! USE m_time_screen                                                                                              !NOADJ
   USE m_numeric
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in   )  ::  mesh
   type( unk ), intent(in   )  ::  dof0
   type( unk ), intent(inout)  ::  dof

   real(rp), intent(out)  ::  cost

   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  sub_nt

   integer(ip)  ::  index_display

   integer(ip)  ::  num_points_display

   !===================================================================================================================!
   !  Model Loop Time Initialization
   !===================================================================================================================!

   dof%h  =  dof0%h
   dof%u  =  dof0%u
   dof%v  =  dof0%v

   tc  =  tc0
   nt  =  nt0

   ! Set serie result for plotting data
! ======================================
   num_points_display = int(ts/dtw) + 1
   if ( (num_points_display-1)*dtw < ts ) then
	num_points_display = num_points_display + 1
   endif
 
   allocate(Series%dof_(num_points_display))

   Series%dof_(1) = dof
   Series%dof_(1)%t_display = 0._rp
! ========================================

   if ( use_obs == 1 ) then
      innovation (:)%ind_t    =  1_ip
      innovW(:)%ind_t  =  1_ip
   end if

   cost  =  0._rp

   !===================================================================================================================!
   !  Writing Initial Condition Output File
   !===================================================================================================================!

   !call write_results( dof0 , mesh )                                                                              !NOADJ

   !===================================================================================================================!
   !  Initializing post treatment variables
   !===================================================================================================================!

   call sw_pre_treatment( dof0 , mesh )                                                                           !NOADJ

   !===================================================================================================================!
   !  SW Model Loop Time
   !===================================================================================================================!

   end_time_loop = .false.

   index_display  = 1

   do while( .not. end_time_loop)

      call sub_run_model
   
   end do


   !===================================================================================================================!
   !  Cost Function Calculation using Innovation Vector
   !===================================================================================================================!
   call calc_cost_function( cost , mesh )


CONTAINS


   SUBROUTINE sub_run_model

      implicit none

      sub_nt = 0

      do while ( .not. end_time_loop .and. sub_nt < max_nt_for_adjoint )

         !=============================================================================================================!
         !  Boundary Conditions
         !=============================================================================================================!

         call set_bc( dof , mesh )

         !=============================================================================================================!
         !  Time Control
         !=============================================================================================================!

         call advance_time( dof , mesh )

         !=============================================================================================================!
         !  Reinitialization of post variables
         !=============================================================================================================!

         bc%sum_mass_flux(:) = 0._rp

         !=============================================================================================================!
         !  Time Stepping Performing ( Euler, RK or IMEX )
         !=============================================================================================================!

         select case( temp_scheme )

            case( 'euler' )

               select case( spatial_scheme  )

                  case( 'first' )

                     call euler_time_step_first( dof , mesh )

                  case( 'first_b1' )

                     call euler_time_step_first_b1( dof , mesh )
                                                                                                                 !<NOADJ
                  case( 'first_b2e' )

   !                  call euler_time_step_first_b2e( dof , mesh )

                  case( 'first_b2i' )

    !                 call euler_time_step_first_b2i( dof , mesh )

                  case( 'first_imp' )

      !               call euler_time_step_imp_first( dof , mesh )

                  case( 'muscl' )

           !          call euler_time_step_muscl( dof , mesh )

                  case( 'muscl_b' )

           !          call euler_time_step_muscl( dof , mesh )

                  case( 'muscl_b1' )

       !              call euler_time_step_muscl_b1( dof , mesh )

                  case( 'muscl_b1_b' )

         !            call euler_time_step_muscl_b1( dof , mesh )

                  case( 'muscl_b2' )

           !          call euler_time_step_muscl_b2( dof , mesh )

                  case( 'muscl_b3' )

   !                  call euler_time_step_muscl_b3( dof , mesh )

                  case( 'muscl_imp' )

      !               call euler_time_step_imp_muscl( dof , mesh )                                                !>NOADJ

      !            case default

                 !    call Stopping_Program_Sub( 'Unknow spatial scheme' )

               end select
                                                                                                                 !<NOADJ
            case( 'rk2' )

           !    call rk_time_step( dof , mesh )

            case( 'jin' )

         !      call jin_time_step( dof , mesh )

            case( 'imex' )

         !      call imex_time_step( dof , mesh )

            case( 'imex2' )

         !      call imex2_time_step( dof , mesh )

            case( 'ck' )

            !   call ck_time_step( dof , mesh )                                                                   !>NOADJ

       !     case default

           !    call Stopping_Program_Sub( 'Unknow temporal scheme' )

         end select
         !=============================================================================================================!
         !  Post-processing
         !=============================================================================================================!

         call sw_post_treatment( dof , mesh )

         !=============================================================================================================!
         !  Writing Output Result File
         !=============================================================================================================!

        ! call write_results( dof , mesh )                                                                         !NOADJ

         !=============================================================================================================!
         !  Filling Innovation Vector
         !=============================================================================================================!

         if ( use_obs == 1 ) then

           ! call calc_innovation( dof,mesh )
           ! call calc_innovW( dof,mesh )
         else

            call update_cost_function( dof , cost )

         endif

         sub_nt = sub_nt + 1

         if ( (abs(tc - index_display*dtw) <= dt .and. tc > index_display*dtw) .or. tc == ts ) then
    	    
	    index_display = index_display + 1

            Series%dof_(index_display) = dof

            Series%dof_(index_display)%t_display = tc

         endif

      end do

   END SUBROUTINE sub_run_model


END SUBROUTINE run_model
