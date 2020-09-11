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
!  Module m_model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_model

   USE m_common
   USE m_linear_algebra
   USE m_mesh
  ! USE m_mpi
  ! USe m_time_screen                                                                                              !NOADJ

   implicit none

   integer(ip)  ::  sw_nb = 3    ! Usefull to create arrays dimensioned with the number of unknows

   !===================================================================================================================!
   !  Discrete Model Unknows Structure
   !===================================================================================================================!

   TYPE unk

      real(rp), dimension(:), allocatable  ::  h
      real(rp), dimension(:), allocatable  ::  u
      real(rp), dimension(:), allocatable  ::  v
      real(rp)     			   ::  t_display

      type( vec2d ), dimension(:), allocatable  ::  grad_h
      type( vec2d ), dimension(:), allocatable  ::  grad_u
      type( vec2d ), dimension(:), allocatable  ::  grad_v
      type( vec2d ), dimension(:), allocatable  ::  grad_z

   END TYPE unk

   TYPE serie

      type(unk), dimension(:), allocatable  ::  dof_

      integer(ip)			    ::  size_series

   END TYPE serie

   type(serie)  ::  Series
   !===================================================================================================================!
   !  Discrete Variables specific to Model
   !===================================================================================================================!

   real(rp), dimension(:), allocatable  ::  bathy_node         ! Bathymetry at mesh nodes
   real(rp), dimension(:), allocatable  ::  bathy_cell         ! Bathymetry at mesh cells gravity center
   real(rp), dimension(:), allocatable  ::  manning            ! Manning coefficient for mesh cells
   integer(ip)  ::  nland                                      ! Total number of land associated to Manning

   integer(ip), dimension(:), allocatable  ::  land            ! Cells land number associated to Manning

   type( vec2d ), dimension(:), allocatable  ::  grad_z        ! Cell Gradient of bathy_cell
   type( vec2d ), dimension(:), allocatable  ::  grad_z2       ! Cell Gradient of bathy_cell^2
   type( vec2d ), dimension(:), allocatable  ::  z_eq          ! Equivalent Bathymetry

   real(rp)  ::  mass_cut
   integer(ip)  ::  manning_data_glob
   !===================================================================================================================!
   !  Boundary Condition Structures
   !===================================================================================================================!

   integer(ip)  ::  feedback_inflow
   real(rp)  ::  coef_feedback

   TYPE hydrograph

      integer(ip)  ::  group

      real(rp), dimension(:), allocatable  ::  t , q

   END TYPE

   TYPE ratcurve

      integer(ip)  ::  group

      real(rp), dimension(:), allocatable  ::  h , q

      real(rp)  ::  z_rat_ref , zout , c1 , c2 , pow(2)

   END TYPE

   TYPE bcs

      integer(ip)  ::  nb , nb_in , nb_out

      character(len=lchar), dimension(:,:), allocatable  ::  typ

      integer(ip), dimension(:), allocatable  ::   grpf

      real(rp), dimension(:), allocatable  ::  inflow
      real(rp), dimension(:), allocatable  ::  outflow

      type(hydrograph), dimension(:), allocatable  ::  hyd

      type(ratcurve), dimension(:), allocatable  ::  rat

      real(rp), dimension(:), allocatable  ::  sum_mass_flux

   END TYPE bcs

   type( bcs )  ::  bc

   !===================================================================================================================!
   !  Recording Structures
   !===================================================================================================================!

   TYPE station_obs

      type( point_in_mesh ), dimension(:), allocatable  ::  pt
      
      real(rp)  :: weight    ! weight of observations

      real(rp)  :: length    ! Length of river
      
      real(rp)  :: dt_offset ! Time of first observation

      real(rp)  :: dt        ! Frequency of observation ( satellite time repetitiveness) 

      real(rp),dimension(:), allocatable :: dt_obs   ! Array with observation time

      integer(ip) :: ind_t   ! Index observation time

      integer(ip) :: nb_dt   ! Number of observation time
      
      real(rp), dimension(:), allocatable  ::  t , h , u , v , q, w

   END TYPE station_obs

   TYPE section_obs

      type( point_in_mesh ), dimension(:), allocatable  ::  pt

      real(rp)  ::  dt , dx

      real(rp), dimension(:), allocatable  ::  t , h , u , v , q

      type( vec2d )  ::  normal

   END TYPE section_obs

   !===================================================================================================================!
   !  Recording Variables
   !===================================================================================================================!

   type( station_obs ), dimension(:), allocatable  ::  station
   type( section_obs ), dimension(:), allocatable  ::  section

   !===================================================================================================================!
   !  Input variables specific to model (in addition to m_common)
   !===================================================================================================================!

   real(rp)     ::  g                                 ! Gravity constant
   real(rp)     ::  heps                              ! Cut-off of water depth to stabilize numerical scheme
   integer(ip)  ::  friction                          ! Activation of a Friction Law in Model
   integer(ip)  ::  c_manning                         !
   integer(ip)  ::  c_bathy                           ! Variables in control vector ( X )
   integer(ip)  ::  c_ic                              !
   integer(ip)  ::  c_hydrograph                      ! derivated by adjoint model ( grad(X) )
   integer(ip)  ::  c_ratcurve                        !

   real(rp)     ::  eps_manning                       !
   real(rp)     ::  eps_bathy                         ! Each variable eps in perturbation control vector
   real(rp)     ::  eps_ic                            !
   real(rp)     ::  eps_hydrograph                    ! used to test validity of the adjoint model
   real(rp)     ::  eps_ratcurve                      !

   real(rp)     ::  regul_manning                     !
   real(rp)     ::  regul_bathy                       !
   real(rp)     ::  regul_ic                          !
   real(rp)     ::  regul_hydrograph                  !
   real(rp)     ::  regul_ratcurve                    !

   integer(ip)  ::  fix_time_step_serie               !

!===================================================================================================================!
   !  Input variables specific to boundary condition
   !===================================================================================================================!


   TYPE boundarycondition

	integer(ip)  ::  nb_in_ , j_hyd

	integer(ip)  ::  nb_out_ , j_rat

	real(rp)     ::  rat_ref

	real(rp), dimension(:), allocatable  ::  t_hyd , q_hyd

	real(rp), dimension(:), allocatable  ::  h_rat , q_rat

   END TYPE boundarycondition

   type(boundarycondition) :: bc_

!===================================================================================================================!
   !  Variables in input.txt
   !===================================================================================================================!

   TYPE inputdata 

	character(len=lchar) :: mesh_type_

	character(len=lchar) :: mesh_name_

	real(rp)	     :: lx_

	real(rp)	     :: ly_

	integer(ip) 	     :: nx_

	integer(ip) 	     :: ny_

	character(len=lchar) :: bc_n_

	character(len=lchar) :: bc_s_

	character(len=lchar) :: bc_w_
 
	character(len=lchar) :: bc_e_

	real(rp)	     :: ts_

	real(rp)	     :: dtw_

	real(rp)	     :: dtp_

	real(rp)	     :: dta_

	character(len=lchar) :: temp_scheme_

	character(len=lchar) :: spatial_scheme_

	integer(ip)	     :: adapt_dt_

	real(rp)	     :: cfl_

	real(rp)	     :: heps_

	integer(ip)	     :: friction_

	real(rp)	     :: g_

	integer(ip)	     :: w_tecplot_

	integer(ip)	     :: w_obs_

	integer(ip)	     :: use_obs_

	real(rp)	     :: eps_min_

	integer(ip)	     :: c_manning_

	real(rp)	     :: eps_manning_

	real(rp)	     :: regul_manning_

	integer(ip)	     :: c_bathy_

	real(rp)	     :: eps_bathy_

	real(rp)	     :: regul_bathy_

	integer(ip)	     :: c_hydrograph_

	real(rp) 	     :: eps_hydrograph_

	real(rp)	     :: regul_hydrograph_

   END TYPE inputdata


   !===================================================================================================================!
   !  Input variables namelist (m_common + model specific upper ones)
   !===================================================================================================================!

   namelist/list_input/ &
      mesh_type, &
      mesh_name, &
      lx, &
      ly, &
      nx, &
      ny, &
      bc_N, &
      bc_S, &
      bc_W, &
      bc_E, &
      ts, &
      dt, &
      dtw, &
      dtp, &
      dta, &
      cfl, &
      adapt_dt, &

      w_vtk, &
      w_tecplot, &
      w_gnuplot, &
      w_exact, &
      w_norm, &
      w_obs, &

      use_obs,&

      spatial_scheme, &
      temp_scheme, &

      max_nt_for_direct , &
      max_nt_for_adjoint, &

      restart_min, &
      eps_min, &

      g, &
      heps, &
      friction, &
      feedback_inflow, &
      coef_feedback, &
      c_manning, &
      c_bathy, &
      c_ic, &
      c_hydrograph, &
      c_ratcurve, &

      eps_manning, &
      eps_bathy, &
      eps_ic, &
      eps_hydrograph, &
      eps_ratcurve, &

      regul_manning, &
      regul_bathy, &
      regul_ic, &
      regul_hydrograph, &
      regul_ratcurve

CONTAINS

!*****************************************************************************************************
! ******************************************************************************************************
! Set input data
!*****************************************************************************************************
! ******************************************************************************************************

   subroutine set_inputdata(inData)

        implicit none

        type(inputdata), intent(in) :: inData

	mesh_type= inData%mesh_type_

	mesh_name= inData%mesh_name_

	lx= inData%lx_

	ly= inData%ly_

	nx= inData%nx_

	ny= inData%ny_

	bc_N= inData%bc_n_

	bc_S= inData%bc_s_

	bc_W= inData%bc_w_

	bc_E= inData%bc_e_

	ts= inData%ts_

	dtw= inData%dtw_

	dtp= inData%dtp_

	dta= inData%dta_

	temp_scheme= inData%temp_scheme_

	spatial_scheme= inData%spatial_scheme_

	adapt_dt= inData%adapt_dt_

	cfl= inData%cfl_

	heps= inData%heps_

	friction= inData%friction_

	g= inData%g_

	w_tecplot= inData%w_tecplot_

	w_obs= inData%w_obs_

	use_obs= inData%use_obs_

	eps_min= inData%eps_min_

	c_manning= inData%c_manning_

	eps_manning= inData%eps_manning_

	regul_manning= inData%regul_manning_

	c_bathy= inData%c_bathy_

	eps_bathy= inData%eps_bathy_

	regul_bathy= inData%regul_bathy_

	c_hydrograph= inData%c_hydrograph_ 

	eps_hydrograph= inData%eps_hydrograph_

	regul_hydrograph= inData%regul_hydrograph_

    end subroutine

    subroutine inputdata_initialise(inData) 

	implicit none

	type(inputdata), intent(out) :: inData

   end subroutine

!*****************************************************************************************************
! ******************************************************************************************************
! Set boundary condition
!*****************************************************************************************************
! ******************************************************************************************************

   subroutine set_boundarycondition( t1 , q1 , h2 , q2 )

	implicit none

	real(rp), dimension(:), intent(in)  ::  t1 , q1

	real(rp), dimension(:), intent(in)  ::  h2 , q2

	bc_%nb_in_ = 1

	bc_%j_hyd = 241

	bc_%nb_out_ = 1
	
	bc_%j_rat = 101 

	bc_%rat_ref = 0._rp

	bc_%t_hyd = t1

	bc_%q_hyd = q1

	bc_%h_rat = h2

	bc_%q_rat = q2

   end subroutine

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Default values for Input variables namelist
!
!**********************************************************************************************************************!
!**********************************************************************************************************************! 

   SUBROUTINE Default_values

      cfl       =  0.8_rp
      adapt_dt  =  1_ip

      dtw       =  2000.
      dtp       =  60.
      dt	=  5.

      w_tecplot =  1_ip
      w_vtk     =  0_ip
      w_gnuplot =  0_ip

      w_exact   =  0_ip
      w_norm    =  0_ip
      w_obs     =  1_ip

      use_obs   =  0_ip

      spatial_scheme  =  'first_b1'
      temp_scheme     =  'euler'

      max_nt_for_direct   =  100000000_ip
      max_nt_for_adjoint  =  2500_ip

      g         =  9.81_rp
      heps      =  0.001
      friction  =  1_ip

      c_manning     =  0_ip
      c_bathy       =  0_ip
      c_ic          =  0_ip
      c_hydrograph  =  1_ip
      c_ratcurve    =  0_ip

      eps_manning     =  0.2_rp
      eps_bathy       =  0.01_rp
      eps_ic          =  0.1_rp
      eps_hydrograph  =  0.2_rp
      eps_ratcurve    =  0.1_rp

      regul_manning     =  0._rp
      regul_bathy       =  0._rp
      regul_ic          =  0._rp
      regul_hydrograph  =  0._rp
      regul_ratcurve    =  0._rp

      inquire( iolength = length_real ) tc

      tc0  =  0._rp
      nt0  =  0_ip

      feedback_inflow   =	1_ip
      coef_feedback     =	0.1_rp
      verbose  =  0_ip

      restart_min  =  0_ip
      eps_min      =  1.d-4

      fix_time_step_serie  =  0_ip

      is_file_open(:)  =  ''
      file_open_counter  =  0

   END SUBROUTINE Default_values

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   subroutine unk_initialise(dof,mesh)
      implicit none
      type(msh), intent(in)  ::  mesh
      type(unk), intent(out)  ::  dof
      allocate(dof%h(mesh%nc + mesh%ncb))
      allocate(dof%u( mesh%nc + mesh%ncb))
      allocate(dof%v(mesh%nc + mesh%ncb))
      allocate(dof%grad_h(mesh%nc + mesh%ncb))
      allocate(dof%grad_u(mesh%nc + mesh%ncb))
      allocate(dof%grad_v(mesh%nc + mesh%ncb))
      allocate(dof%grad_z(mesh%nc + mesh%ncb))

      dof%h(:)  =  0._rp
      dof%u(:)  =  0._rp
      dof%v(:)  =  0._rp

      dof%grad_h(:)%x  =  0._rp
      dof%grad_h(:)%y  =  0._rp

      dof%grad_u(:)%x  =  0._rp
      dof%grad_u(:)%y  =  0._rp

      dof%grad_v(:)%x  =  0._rp
      dof%grad_v(:)%y  =  0._rp

      dof%grad_z(:)%x  =  0._rp
      dof%grad_z(:)%y  =  0._rp   
   end subroutine


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   SUBROUTINE unk_finalise(dof)
      implicit none
      type(unk), intent(inout)  ::  dof
      if (allocated(dof%h)) deallocate(dof%h)
      if (allocated(dof%u)) deallocate(dof%u)
      if (allocated(dof%v)) deallocate(dof%v)
      if (allocated(dof%grad_h)) deallocate(dof%grad_h)
      if (allocated(dof%grad_u)) deallocate(dof%grad_u)
      if (allocated(dof%grad_v)) deallocate(dof%grad_v)
      if (allocated(dof%grad_z)) deallocate(dof%grad_z)
   END SUBROUTINE 


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE dealloc_model

      implicit none

      if ( allocated( bathy_node ) ) deallocate( bathy_node )
      if ( allocated( bathy_cell ) ) deallocate( bathy_cell )
      if ( allocated( land       ) ) deallocate( land       )
      if ( allocated( manning    ) ) deallocate( manning    )
   END SUBROUTINE dealloc_model



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calling MPI and communicating dof to fill ghost cells
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


!   SUBROUTINE com_dof( dof , mesh )

 !     implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

  !    type( msh ), intent(in   )  ::  mesh
   !   type( unk ), intent(inout)  ::  dof

      !================================================================================================================!
      !
      !================================================================================================================!

    !  #if defined USE_MPI || USE_MPI_ADJ

     !    call Time_Init_Part(80)                                                                                  !NOADJ

     !    call com_var_r( dof%h(:) , mesh )
      !   call com_var_r( dof%u(:) , mesh )
       !  call com_var_r( dof%v(:) , mesh )

        ! call Time_End_Part(80)                                                                                   !NOADJ

     ! #endif

  ! END SUBROUTINE com_dof


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation or Reallocate Stations Type
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

  ! SUBROUTINE alloc_or_realloc_station( station_inout , new )

  !    implicit none

  !    type( station_obs ), dimension(:), allocatable, intent(inout)  ::  station_inout

  !    integer(ip), intent(in)  ::  new

 !     integer(ip)  ::  old , iobs , pt

  !    type( station_obs ), dimension(:), allocatable  ::  station_tmp

 !     intrinsic move_alloc

   !   if ( .not. allocated( station_inout ) ) then

   !      allocate( station_inout( new ) )

      !   return

   !   end if

   !   old = size( station_inout )

    !  if     ( new == old ) then

   !      return

   !   else if ( new > old ) then

   !      allocate( station_tmp( new ) )

   !      do iobs = 1,old

     !       station_tmp( iobs )%dt      =  station_inout( iobs )%dt
     !       station_tmp( iobs )%weight  =  station_inout( iobs )%weight

      !      allocate( station_tmp( iobs )%pt( size(station_inout( iobs )%pt ) ) )

       !     do pt = 1,size( station_inout( iobs )%pt )

       !        station_tmp( iobs )%pt( pt )%cell   =  station_inout( iobs )%pt( pt )%cell
        !       station_tmp( iobs )%pt( pt )%coord  =  station_inout( iobs )%pt( pt )%coord

        !    end do

       !  end do

     !    call move_alloc( station_tmp , station_inout )

   !   else

    !     call Stopping_Program_Sub( 'Wrong Station Dimension for Allocation' )

   !   end if

  ! END SUBROUTINE alloc_or_realloc_station

END MODULE m_model
