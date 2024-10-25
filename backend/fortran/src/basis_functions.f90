#include "numbat_decl.h"

module class_BasisFunctions

   use numbatmod
   use alloc

   type, public :: BasisFunctions


      ! phi_P1_ref = values of the P2 polynomials (1-6) at a given point in the reference triangle.
      ! gradt_P1_ref = gradients of the P2s at that point
      ! gradt_P1_act = gradients of the P2s in the actual triangle
      double precision phi_P1_ref(P1_NODES_PER_EL)
      double precision gradt_P1_ref(2, P1_NODES_PER_EL)
      double precision gradt_P1_act(2, P1_NODES_PER_EL)


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

      integer(8) phi_vec_map(4,3,N_DDL_T)


   contains

   procedure :: make_phi_vector_map => BasisFunctions_make_phi_vector_map
   procedure :: evaluate_at_position => BasisFunctions_evaluate_at_position
   end type BasisFunctions

   contains

#include "basis_functions_impl.f90"

end module
