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

      double precision phi_P2_ref(P2_NODES_PER_EL)
      double precision gradt_P2_ref(2, P2_NODES_PER_EL)
      double precision gradt_P2_act(2, P2_NODES_PER_EL)


      double precision phi_P3_ref(P3_NODES_PER_EL)
      double precision gradt_P3_ref(2, P3_NODES_PER_EL)
      double precision gradt_P3_act(2, P3_NODES_PER_EL)


      double precision mat_B(2,2)    !transform ref to act triangle
      double precision mat_T(2,2)    ! inverse transform T=B^-1
      double precision det   ! transformation determinant

   contains

      procedure :: setup_reference_quadratures => QuadIntegrator_setup_reference_quadratures

      procedure :: build_transforms_at => QuadIntegrator_build_transforms_at
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

      integer(8) n_msh_el
      integer(8) n_msh_pts
      integer(8), dimension(:,:), allocatable :: elnd_to_mesh
      double precision, dimension(:,:), allocatable :: v_nd_xy

   contains
      procedure :: init_from_py => PyFrontEnd_init_from_py
      procedure :: nodes_at_el => PyFrontEnd_nodes_at_el

   end type PyFrontEnd

contains

#include "quad_integrator_impl.f90"

#include "analytic_integrator_impl.f90"


end module class_TriangleIntegrators
