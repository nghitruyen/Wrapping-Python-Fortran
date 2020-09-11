!======================================================================================================================!
!
!                    DassFlow1D Version 2.0
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
!    - 1D Shallow Module (1D Shallow Water Model, Finite Volume/Finite Difference Methods), i.e. the present code.
!    - 2D Shallow Module (2D Shallow Water Model, Finite Volume Method)
!    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current main developers or scientific contributers are:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               K. Larnier (C.S Communication and Systems & INSA Toulouse).
!               J. Monnier (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA Toulouse & Mathematics Institute of Toulouse IMT).
!               P.-A. Garambois (INSA Strasbourg & ICUBE).
!               L. Pujol (CNES & INSA Strasbourg & ICUBE).
!  and former other developers (P. Brisset, R. Madec, M. Honnorat and J. Marin).
!
!  Scientific Contact : jerome.monnier@insa-toulouse.fr
!  Technical  Contact : kevin.larnier@c-s.fr
!                       frederic.couderc@math.univ-toulouse.fr
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
!> @file m_sw_mono.f90
!> @brief This file includes the m_model module.

!> @brief Module m_sw_mono.
module m_sw_mono
    use m_common
    use m_linear_algebra
    use m_mesh, only: Mesh
    use m_obs
!     use m_time_screen                                                                                              !NOADJ
    implicit none
#ifndef CPP_ADJ
    include "mpif.h"
    include 'dmumps_struc.h'
#endif
   
   !===================================================================================================================!
   !  Discrete Model Unknows Structure
   !===================================================================================================================!

    !> @brief Structure of unknowns.
    !> @details This structure includes all the unknows of models (A, Q, h, qlat) and the terms for compuing the missfit
    type Unknowns

        !> Flow areas
        real(rp), dimension(:), allocatable :: A
        !> Discharges
        real(rp), dimension(:), allocatable :: Q
        !> Depths
        real(rp), dimension(:), allocatable :: h
        !> Lateral inflows
        real(rp), dimension(:,:), allocatable :: qlat
        !> Friction term values
        real(rp), dimension(:), allocatable :: Sf
        !> Gravity term values
        real(rp), dimension(:), allocatable :: Sg

    end type

    !> Structure of implicit matrix
    type ImplicitMatrix

        !> Segments offsets
        integer(ip), dimension(:), allocatable :: seg_offsets
        !> GA values
        real(rp), dimension(:), allocatable :: GA
        !> GB values
        real(rp), dimension(:), allocatable :: GB
        !> GC values
        real(rp), dimension(:), allocatable :: GC
        !> GD values
        real(rp), dimension(:), allocatable :: GD
        !> GE values
        real(rp), dimension(:), allocatable :: GE
        !> GF values
        real(rp), dimension(:), allocatable :: GF
        !> CR values
        real(rp), dimension(:), allocatable :: CR
        !> CS values
        real(rp), dimension(:), allocatable :: CS
        !> CT values
        real(rp), dimension(:), allocatable :: CT
        
        ! TA values (condensed relations x1 = TA.xn + TB)
        real(rp), dimension(:), allocatable :: TA1
        real(rp), dimension(:), allocatable :: TA2
        real(rp), dimension(:), allocatable :: TA3
        real(rp), dimension(:), allocatable :: TA4
        ! TB values (condensed relations x1 = TA.xn + TB)
        real(rp), dimension(:), allocatable :: TB1
        real(rp), dimension(:), allocatable :: TB2

#ifndef CPP_ADJ     
        ! MUMPS System (condensed system)
        type(DMUMPS_STRUC) :: sys
#endif
        real(rp), dimension(:), allocatable :: ANZ
        real(rp), dimension(:), allocatable :: RHS

    end type

    !> Structure of timeseries.
    type Timeseries
        real(rp), dimension(:), allocatable  ::  t  !< Time array.
        real(rp), dimension(:), allocatable  ::  y  !< Flow array.
    end type

    !> Structure of boundary condition
    type BoundaryCondition

        !> ID
        character(len=16) :: id
        !> Timeseries
        type(Timeseries) :: ts

    end type

    !> Structure of inflow condition
    type InflowCondition

        !> Index of segment
        integer(ip) :: iseg
        !> Index of edge (between cross-sections [ie, ie+1])
        integer(ip) :: ie
        !> Timeseries
        type(Timeseries) :: ts

    end type
    
    type Results
    
        real(rp), dimension(:), allocatable :: t
        real(rp), dimension(:, :), allocatable :: q
        real(rp), dimension(:, :), allocatable :: h
        real(rp), dimension(:, :), allocatable :: a
    
    end type
        
    type Model

        !> Flag to enable discharge in estimations
        logical :: discharge_estimation
        !> Status code
        integer(ip) :: status
        !> Warning counters (1=number of heps correction encountered, 2=TBD, 3=TBD, 4=TBD, 5=TBD)
        integer(ip), dimension(5) :: warning_counters
        !> Current time
        real(rp) :: tc
        !> End time
        real(rp) :: te
        !> Start time
        real(rp) :: ts
        !> Timestep
        real(rp) :: dt
        !> Timestep for printing informations on standard output
        real(rp) :: dtout
        !> Froude threshold for LPI method
        real(rp) :: frLPI
        !> Gravity acceleration (m2/s)
        real(rp) :: gravity
        !> Depth threshold
        real(rp) :: heps
        !> Discharge threshold
        real(rp) :: qeps
        !> Parameter for LPI method
        real(rp) :: mLPI
        !> Implicit coefficient for Preissmann scheme
        real(rp) :: theta_preissmann
        !> Regularization coefficient
!         real(rp) :: alpha_reg
        real(rp) :: gamma_reg
        !> Observation cost
        real(rp) :: cost_obs
        !> Regularization cost
        real(rp) :: cost_reg
        !> Computational scheme
        character(len=32) :: scheme
        !>scale for the implicit diffusive wave
        character(len=16) :: scale_idw
        !> Bounday conditions
        type(BoundaryCondition), dimension(:), allocatable :: bc
        !> Inflow conditions
        type(InflowCondition), dimension(:), allocatable :: ic
        !> Mesh
        type(Mesh), pointer :: msh
        !> large grid
        type(Mesh), pointer :: large_grid
        !> Unknowns
        type(Unknowns) :: dof
        !> Implicit arrays
        type(ImplicitMatrix) :: imp
        !> Results
        type(Results) :: res
      
    end type


#ifndef CPP_ADJ
    interface
        subroutine generate_observations(mdl, obs)
            import Model
            import Observations
            implicit none
            type(Model), intent(inout) :: mdl
            type(Observations), intent(inout), optional :: obs
        end subroutine
        subroutine standard_step(mdl, eps, itermax)
            import ip
            import Model
            import rp
            implicit none
            type(Model), intent(inout) :: mdl
            real(rp), intent(in) :: eps
            integer(ip), intent(in) :: itermax
        end subroutine
        subroutine preissmann_timestep(mdl, msh, imp, dof, status)
            import ip
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            integer(ip) :: status
        end subroutine
	subroutine implicit_diffusive_wave(mdl, msh, large_grid, imp, dof, status)
            import ip
            import ImplicitMatrix
            import Mesh
            import Model
            import Unknowns
            implicit none
            type(Model), intent(in) :: mdl
            type(Mesh), intent(inout) :: msh
            type(Mesh), intent(inout) :: large_grid
            type(ImplicitMatrix), intent(inout) :: imp
            type(Unknowns), intent(inout) :: dof
            integer(ip) :: status
        end subroutine
        subroutine time_loop(mdl, obs)
            import Model
            import Observations
            implicit none
            type(Model), intent(inout) :: mdl
            type(Observations), intent(inout), optional :: obs
        end subroutine
    end interface
#endif

      
#ifndef CPP_ADJ
    contains
    
    
    subroutine set_timeseries(bc, t, y)
        implicit none
        type(BoundaryCondition), intent(inout) :: bc
        real(rp), dimension(:), intent(in) :: t
        real(rp), dimension(:), intent(in) :: y

        if (size(t) /= size(y)) then
            call f90wrap_abort("t and y must be of same size")
        end if
        
        if (allocated(bc%ts%t)) deallocate(bc%ts%t)
        if (allocated(bc%ts%y)) deallocate(bc%ts%y)
        
        allocate(bc%ts%t(size(t)))
        allocate(bc%ts%y(size(t)))
        bc%ts%t(:) = t
        bc%ts%y(:) = y
        
    end subroutine
    
    
    subroutine implicitmatrix_initialise(imp, mdl, msh)
        implicit none
        type(Model), intent(in) :: mdl
        type(Mesh), intent(in) :: msh
        type(ImplicitMatrix), intent(inout) :: imp
        
        ! Index of segment
        integer(ip) :: iseg
        ! Number of rows
        integer(ip) :: nrows
        ! Number of boundary conditions
        integer(ip) :: nbc
        
        ! Allocate offsets array
        allocate(imp%seg_offsets(msh%nseg))
        
        ! Compute offsets and total number of rows
        nrows = 0
        do iseg = 1, msh%nseg
        
            imp%seg_offsets(iseg) = nrows
            nrows = nrows + msh%seg(iseg)%last_cs - msh%seg(iseg)%first_cs + 1

        end do
        
        ! Allocate coefficients arrays
        allocate(imp%GA(nrows))
        allocate(imp%GB(nrows))
        allocate(imp%GC(nrows))
        allocate(imp%GD(nrows))
        allocate(imp%GE(nrows))
        allocate(imp%GF(nrows))
        allocate(imp%CR(nrows))
        allocate(imp%CS(nrows))
        allocate(imp%CT(nrows))
        
        ! Allocate coefficients arrays for condensed relations
        allocate(imp%TA1(msh%nseg))
        allocate(imp%TA2(msh%nseg))
        allocate(imp%TA3(msh%nseg))
        allocate(imp%TA4(msh%nseg))
        allocate(imp%TB1(msh%nseg))
        allocate(imp%TB2(msh%nseg))
        
        call init_implicit_reduced_system(mdl, msh, imp)
        
    end subroutine
    
    
    subroutine implicitmatrix_finalise(imp)
        implicit none
        type(ImplicitMatrix), intent(inout) :: imp
        
        ! Deallocate arrays
        if (allocated(imp%seg_offsets)) deallocate(imp%seg_offsets)
        if (allocated(imp%GA)) deallocate(imp%GA)
        if (allocated(imp%GB)) deallocate(imp%GB)
        if (allocated(imp%GC)) deallocate(imp%GC)
        if (allocated(imp%GD)) deallocate(imp%GD)
        if (allocated(imp%GE)) deallocate(imp%GE)
        if (allocated(imp%GF)) deallocate(imp%GF)
        if (allocated(imp%CR)) deallocate(imp%CR)
        if (allocated(imp%CS)) deallocate(imp%CS)
        if (allocated(imp%CT)) deallocate(imp%CT)
        if (allocated(imp%TA1)) deallocate(imp%TA1)
        if (allocated(imp%TA2)) deallocate(imp%TA2)
        if (allocated(imp%TA3)) deallocate(imp%TA3)
        if (allocated(imp%TA4)) deallocate(imp%TA4)
        if (allocated(imp%TB1)) deallocate(imp%TB1)
        if (allocated(imp%TB2)) deallocate(imp%TB2)
        
    end subroutine
    
    
    subroutine model_initialise(mdl, msh, scheme)
        implicit none
        type(Mesh), intent(in), target :: msh
        type(Model), intent(out) :: mdl
        character(len=*), intent(in), optional :: scheme
        
        ! Index of boundary condition
        integer(ip) :: ibc
        ! Index of segment
        integer(ip) :: iseg
        ! Number of boundary conditions
        integer(ip) :: nbc

        ! Set mesh pointer
        mdl%msh => msh
        ! Set large grid default pointer
        mdl%large_grid => NULL()       
        
        ! Compute number of boundary conditions
        nbc = 0
        do iseg = 1, msh%nseg
        
            if (msh%seg(iseg)%us_bc > 0) then
                nbc = max(nbc, msh%seg(iseg)%us_bc)
            end if
            if (msh%seg(iseg)%ds_bc > 0) then
                nbc = max(nbc, msh%seg(iseg)%ds_bc)
            end if

        end do
        
        ! Allocate and initialise boundary conditions
        allocate(mdl%bc(nbc))
        do ibc = 1, nbc
            mdl%bc(ibc)%id = "undefined"
        end do
        
        ! Set default values
        mdl%discharge_estimation = .false.
        mdl%tc = 0.0
        mdl%ts = 0.0
        mdl%te = 1.0
        mdl%dt = 1.0
        mdl%dtout = 1.0
        mdl%heps = 0.01
        mdl%qeps = 0.001
        mdl%frLPI = 0.7
        mdl%mLPI = 10.0
        mdl%gravity = 9.81
        mdl%theta_preissmann = 0.7
!         mdl%alpha_reg = 0.0
        mdl%gamma_reg = -1.0
        mdl%scale_idw = "low"
        
        ! Allocate unknowns
        call unknowns_initialise(mdl%dof, mdl%msh)
        
        mdl%scheme = "undefined"
        
    end subroutine
    
    subroutine set_large_grid(mdl,large_grid) 

	implicit none
        type(Model), intent(inout) :: mdl
        type(Mesh), intent(in), target :: large_grid
        
        mdl%large_grid => large_grid

    end subroutine    

    subroutine add_inflow_condition(mdl, iseg, coords, t, q)
        implicit none
        type(Model) :: mdl
        integer(ip) :: iseg
        real(rp), dimension(2) :: coords
        real(rp), dimension(:) :: t
        real(rp), dimension(:) :: q
        
        ! Iterator
        integer(ip) :: i
        ! Index of cross-section
        integer(ip) :: ics
        ! Index of edge
        integer(ip) :: ie
        ! Index of inflow condition
        integer(ip) :: iic
        ! Number of inflow conditions before adding new inflow condition
        integer(ip) :: size_old
        ! Squares of distances
        real(rp), dimension(2) :: dist2
        ! Minimum square of distance
        real(rp) :: dist2min
        ! Temporary arry for copying existing inflow conditions
        type(InflowCondition), dimension(:), allocatable :: old_ic
        
        ! CHECK-UP
        if (iseg < 0 .or. iseg > size(mdl%msh%seg)) then
            call f90wrap_abort("index of segment out of bounds")
        end if
        if (size(t) /= size(q)) then
            call f90wrap_abort("t and q must be of same size")
        end if

        ! Locate closest cross-section
        ie = mdl%msh%seg(iseg+1)%first_cs
        ics = ie
        dist2min = (coords(1) - mdl%msh%cs(ics)%coord%x)**2 + (coords(2) - mdl%msh%cs(ics)%coord%y)**2
        do ics = mdl%msh%seg(iseg+1)%first_cs+1, mdl%msh%seg(iseg+1)%last_cs
            dist2(1) = (coords(1) - mdl%msh%cs(ics)%coord%x)**2 + (coords(2) - mdl%msh%cs(ics)%coord%y)**2
            if (dist2min > dist2(1)) then
                dist2min = dist2(1)
                ie = ics
            end if
        end do
        
        ! Update index of edge
        if (ie > mdl%msh%seg(iseg+1)%first_cs .and. ie < mdl%msh%seg(iseg+1)%last_cs) then
            dist2(1) = (coords(1) - mdl%msh%cs(ie-1)%coord%x)**2 + (coords(2) - mdl%msh%cs(ie-1)%coord%y)**2
            dist2(2) = (coords(1) - mdl%msh%cs(ie+1)%coord%x)**2 + (coords(2) - mdl%msh%cs(ie+1)%coord%y)**2
            if (dist2(1) < dist2(2)) ie = ie - 1
        end if
        
        
        ! Check that inflow conditions are set in order (increasing iseg, ics)
        if (allocated(mdl%ic)) then
            size_old = size(mdl%ic)
            if (iseg+1 < mdl%ic(size_old)%iseg) then
                call f90wrap_abort("inflow conditions must be added with increasing indices of segments")
            end if
            if (iseg+1 == mdl%ic(size_old)%iseg .and. ie < mdl%ic(size_old)%ie) then
                call f90wrap_abort("inflow conditions must be added with increasing indices of edges in a segment")
            end if
        end if
        
        
        ! Resize array of inflow conditions
        if (allocated(mdl%ic)) then
            size_old = size(mdl%ic)
            allocate(old_ic(size(mdl%ic)))
            do i = 1, size(mdl%ic)
                old_ic(i)%iseg = mdl%ic(i)%iseg
                old_ic(i)%ie = mdl%ic(i)%ie
                allocate(old_ic(i)%ts%t(size(mdl%ic(i)%ts%t)))
                old_ic(i)%ts%t(:) = mdl%ic(i)%ts%t(:)
                allocate(old_ic(i)%ts%y(size(mdl%ic(i)%ts%y)))
                old_ic(i)%ts%y(:) = mdl%ic(i)%ts%y(:)
                deallocate(mdl%ic(i)%ts%t)
                deallocate(mdl%ic(i)%ts%y)
            end do
            deallocate(mdl%ic)
        else
            size_old = 0
        end if
        allocate(mdl%ic(size_old+1))
        if (size_old > 0) then
            do i = 1, size_old
                mdl%ic(i)%iseg = old_ic(i)%iseg
                mdl%ic(i)%ie = old_ic(i)%ie
                allocate(mdl%ic(i)%ts%t(size(old_ic(i)%ts%t)))
                mdl%ic(i)%ts%t(:) = old_ic(i)%ts%t(:)
                allocate(mdl%ic(i)%ts%y(size(old_ic(i)%ts%t)))
                mdl%ic(i)%ts%y(:) = old_ic(i)%ts%y(:)
                deallocate(old_ic(i)%ts%t)
                deallocate(old_ic(i)%ts%y)
            end do
            deallocate(old_ic)
        end if
        
        ! Set new inflow condition
        mdl%ic(size_old+1)%iseg = iseg + 1
        mdl%ic(size_old+1)%ie = ie
        allocate(mdl%ic(size_old+1)%ts%t(size(t)))
        mdl%ic(size_old+1)%ts%t(:) = t(:)
        allocate(mdl%ic(size_old+1)%ts%y(size(q)))
        mdl%ic(size_old+1)%ts%y(:) = q(:) / mdl%msh%cs(ie)%deltademi
        
    end subroutine
    
    
    subroutine set_scheme(mdl, scheme)
        implicit none
        type(Model), intent(inout) :: mdl
        character(len=*), intent(in), optional :: scheme
        
        ! Index of boundary condition
        integer(ip) :: ibc
        ! Index of segment
        integer(ip) :: iseg
        ! Number of boundary conditions
        integer(ip) :: nbc
    
        ! CHECK-UP
        do ibc = 1, size(mdl%bc)
        
            if (mdl%bc(ibc)%id == "undefined") then
                call abort_solver("All boundary conditions must be set before calling 'set_scheme'")
            end if
        
        end do
        ! Allocate implicit arrays
        if (present(scheme)) then
            if (scheme == "preissmann") then
                mdl%scheme = "preissmann"
                call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
                
            end if
	    if (scheme == "implicit_diffusive_wave") then
                mdl%scheme = "implicit_diffusive_wave"
                call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
            end if
            
        else
            mdl%scheme = "implicit_diffusive_wave"
            !mdl%scheme = "preissmann"
            call implicitmatrix_initialise(mdl%imp, mdl, mdl%msh)
        
        end if
        
    end subroutine
    
    
    subroutine model_finalise(mdl)
        implicit none
        type(Model), intent(inout) :: mdl
        
        ! Deallocate boundary conditions
        if (allocated(mdl%bc)) then
            ! TODO finalise bcs
            deallocate(mdl%bc)
        end if
        
    end subroutine
    
    
    subroutine timeseries_initialise(ts, nt, y0)
        implicit none
        type(Timeseries), intent(out) :: ts
        integer(ip), intent(in) :: nt
        real(rp), optional :: y0
        
        if (nt == 0) then
            allocate(ts%y(1))
        else
            allocate(ts%t(nt))
            allocate(ts%y(nt))
        end if
        if (present(y0)) then
            ts%y(:) = y0
        else
            ts%y(:) = 0.0_rp
        end if
        
    end subroutine
    
    
    subroutine timeseries_finalise(ts)
        implicit none
        type(Timeseries), intent(inout) :: ts
        
        ! Deallocate array
        if (allocated(ts%t)) deallocate(ts%t)
        if (allocated(ts%y)) deallocate(ts%y)
        
    end subroutine
    

     !> @brief Initialise the unknowns of the model.
    subroutine unknowns_initialise(dof, msh)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Unknowns
        type(Unknowns), intent(inout) :: dof
        !> Mesh
        type(Mesh), intent(in ) :: msh

        ! Allocate arrays
        allocate(dof%Q(msh%ncs))
        allocate(dof%A(msh%ncs))
        allocate(dof%h(msh%ncs))
        allocate(dof%sg(msh%ncs))
        allocate(dof%sf(msh%ncs))


        ! Set initial values
        dof%Q(:)  =  0._rp
        dof%A(:)  =  0._rp
        dof%h(:)  =  0._rp
        dof%Sg(:)  =  0._rp
        dof%Sf(:) =  0._rp
        
    end subroutine


    subroutine unknowns_finalise(dof)
        implicit none
        !=ARGUMENTS====================================================================================================!
        !> Unknowns
        type(Unknowns), intent(inout) :: dof
        
        if (allocated(dof%Q)) deallocate(dof%Q)
        if (allocated(dof%A)) deallocate(dof%A)
        if (allocated(dof%h)) deallocate(dof%h)
        if (allocated(dof%Sf)) deallocate(dof%Sf)
        if (allocated(dof%Sg)) deallocate(dof%Sg)
        
    end subroutine    
#endif

end module m_sw_mono
