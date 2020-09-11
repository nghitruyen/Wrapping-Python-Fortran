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
!  Module call_run_model : This file is for wrapper part.
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE call_model

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_model
   USE m_linear_solver
   USE m_numeric
   USE m_obs

   implicit none

   TYPE Model

	type(msh), allocatable  ::  mesh

   	type(unk), allocatable  ::  dof0

   	type(unk), allocatable  ::  dof

   	real(rp), allocatable ::  cost

   END TYPE

   TYPE points

	real(rp), dimension(:), allocatable :: x_space

	real(rp), dimension(:), allocatable :: y_space

   END TYPE

CONTAINS

   subroutine Model_initialise(mdl) 

	implicit none

	type(Model), intent(out) :: mdl
	
	allocate(mdl%mesh)
	allocate(mdl%dof0)
	allocate(mdl%dof)
	allocate(mdl%cost)
	mdl%cost = 0._rp

   end subroutine

   subroutine Model_finalise(mdl)

	implicit none

	type(Model), intent(inout) :: mdl

	if (allocated(mdl%mesh)) deallocate(mdl%mesh)
	if (allocated(mdl%dof0)) deallocate(mdl%dof0)
	if (allocated(mdl%dof)) deallocate(mdl%dof)
	if (allocated(mdl%cost)) deallocate(mdl%cost)

   end subroutine

   subroutine run_direct(mdl,t1 , q1 , h2 , q2, inData, res,cost_func,pts)

	implicit none

	type(Model), intent(inout) :: mdl

	type(serie), intent(out) :: res

	real(rp), intent(out) :: cost_func

	type(points), intent(out) :: pts 

	real(rp), dimension(:), intent(in)  ::  t1 , q1

	real(rp), dimension(:), intent(in)  ::  h2 , q2

	type(inputdata), intent(in) :: inData

	call set_boundarycondition( t1 , q1 , h2 , q2 ) ! Set boundary conditions

        call Default_values

	call set_inputdata(inData) ! Set input data

        call reading_args

        call Machine_Number_Limits

	call Init_Linear_Solver(mdl%mesh)
	
	call Mesh_Geometric_Properties(mdl%mesh)

	call Initial(mdl%dof0,mdl%mesh)

	call Init_Schemes(mdl%mesh)

	call run_model(mdl%mesh, mdl%dof0, mdl%dof, mdl%cost)
 
        res%dof_ = Series%dof_

	res%size_series = size(Series%dof_)

	allocate(pts%x_space(mdl%mesh%nc))

	allocate(pts%y_space(mdl%mesh%nc))

	! Set space
        do k = 1,mdl%mesh%nc

             pts%x_space(k) = mdl%mesh%node(k)%coord%x

             pts%y_space(k) = mdl%mesh%node(k)%coord%y

        end do

	cost_func = mdl%cost

   end subroutine
	

END MODULE call_model
