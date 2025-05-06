#include "numbat_decl.h"

module class_TriangleIntegrators

   use numbatmod
   use alloc

   integer(8), parameter :: NQUAD_MAX = 16

   private

   ! A class for 2D numerical integration on the unit triangle

   type, public  :: QuadIntegrator

      integer(8) n_quad  ! number of sample points, usually 16
      double precision x_quad(NQUAD_MAX)
      double precision y_quad(NQUAD_MAX)
      double precision wt_quad(NQUAD_MAX)


      ! phi_P2_ref = values of the P2 polynomials (1-6) at a given point in the reference triangle.
      ! gradt_P2_ref = gradients of the P2s at that point
      ! gradt_P2_act = gradients of the P2s in the actual triangle
      double precision phi_P2_ref(P2_NODES_PER_EL)
      double precision gradt_P2_ref(2, P2_NODES_PER_EL)
      double precision gradt_P2_act(2, P2_NODES_PER_EL)

      ! phi_P3_ref = values of the P3 polynomials (1-10) at a given point in the reference triangle.
      ! gradt_P3_ref = gradients of the P3s at that point
      ! gradt_P3_act = gradients of the P3s in the actual triangle
      double precision phi_P3_ref(P3_NODES_PER_EL)
      double precision gradt_P3_ref(2, P3_NODES_PER_EL)
      double precision gradt_P3_act(2, P3_NODES_PER_EL)

      double precision mat_B(2,2)    ! transform ref to act triangle
      double precision mat_T(2,2)    ! inverse transform T=B^-1
      double precision det           ! transformation determinant
      double precision quadweight    ! transformation determinant

   contains

      procedure :: setup_reference_quadratures => QuadIntegrator_setup_reference_quadratures

      procedure :: build_transforms_at => QuadIntegrator_build_transforms_at

      procedure :: get_current_quadweight => QuadIntegrator_get_current_quadweight

      procedure :: get_quad_point => QuadIntegrator_get_quad_point



   end type QuadIntegrator


   type, public  :: AnalyticIntegrator

      double precision mat_B(2,2)       !transform ref to act triangle
      double precision mat_T(2,2)       ! inverse transform T=B^-1
      double precision mat_T_tr(2,2)    ! inverse transpose transform T_tr= Transpose(B^-1)
      double precision det              ! transformation determinant

   contains

      procedure :: build_transforms_at => AnalyticIntegrator_build_transforms_at


   end type AnalyticIntegrator


   type, public :: PyFrontEnd

      integer(8) n_msh_elts
      integer(8) n_msh_pts
      integer(8), dimension(:,:), allocatable :: m_elnd_to_mshpt
      double precision, dimension(:,:), allocatable :: v_mshpt_xy

      integer(8) cur_elt
      double precision cur_nds_xy(2, P2_NODES_PER_EL)

   contains
      procedure :: init_from_py => PyFrontEnd_init_from_py
      procedure :: init_from_py2 => PyFrontEnd_init_from_py2

      procedure :: nodes_at_el => PyFrontEnd_nodes_at_el
      procedure :: elt_is_curved => PyFrontEnd_elt_is_curved
   end type PyFrontEnd

contains

#include "quad_integrator_impl.f90"

#include "analytic_integrator_impl.f90"


end module class_TriangleIntegrators
