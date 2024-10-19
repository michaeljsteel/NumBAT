#include "numbat_decl.h"

module class_QuadIntegrator

   use numbatmod
   use alloc

   integer(8), parameter :: NQUAD_MAX = 16

   private

   ! A class for 2D numerical integration on the unit triangle

   type, public  :: QuadIntegrator



      integer n_quad  ! number of sample points, usually 16
      double precision x_quad(NQUAD_MAX)
      double precision y_quad(NQUAD_MAX)
      double precision wt_quad(NQUAD_MAX)

      double precision phi_P2_ref(P2_NODES_PER_EL)
      double precision gradt_P2_ref(2, P2_NODES_PER_EL)
      double precision gradt_P2_act(2, P2_NODES_PER_EL)

      double precision mat_B(2,2)    !transform ref to act triangle
      double precision mat_T(2,2)    ! inverse transform T=B^-1
      double precision det   ! transformation determinant

   contains

      procedure :: setup_reference_quadratures => QuadIntegrator_setup_reference_quadratures

      procedure :: build_transforms_at => QuadIntegrator_build_transforms_at
   end type QuadIntegrator

contains

#include "quad_integrator_impl.f90"

end module class_QuadIntegrator
