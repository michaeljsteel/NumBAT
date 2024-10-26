#include "numbat_decl.h"
!  P2 basis function over the unit triangle

!  Quadratic basis function = P2 * Grad P1

!  vector_elt_map(1,j,i) = k : number on data to be stored: if k=3 only one gradient will be used; k=4 => 2 gradients
!  vector_elt_map(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  vector_elt_map(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial phi_n
!  vector_elt_map(4,j,i)     : it will be used only if k=4

!  The basis functions are stored in the 3rd index i as : centre, edges 1(4), 2(5), 3(6)
!                                 in the 2nd index j as 3 functions on each edge (for i=1), and along each edge (for i=2..4)


subroutine BasisFunctions_make_vector_elt_map(this, el_nds)
   class(BasisFunctions) this
   integer(8) el_nds(P2_NODES_PER_EL)     ! global indices of the P2 nodes at current element

   call make_vector_elt_map(el_nds, this%vector_elt_map)
end subroutine

! Evaluates scalar P1, P2 and P3 functions and gradients
! at position t_xy in the reference triangle for an element
! with P2 nodes at el_nds_xy
! Finds the transformation matrices and et for the mapping
! from reference unit  triangle to actual triangle
subroutine  BasisFunctions_evaluate_at_position(this, i_el, t_xy, is_curved, el_nds_xy, nberr)
   use numbatmod

   class(BasisFunctions) this
   double precision t_xy(2), xy_act
   double precision el_nds_xy(2,P2_NODES_PER_EL)
   integer(8) i_el

   logical is_curved
   type(NBError) nberr
   integer debug

   integer(8)  :: errco
   character(len=EMSG_LENGTH)  :: emsg
   errco=0

   debug =0
   !  Elements and gradients for the P1, P2, P3 basis functions
   call phi1_2d_mat (t_xy, this%phi_P1_ref, this%gradt_P1_ref)
   call phi2_2d_mat (t_xy, this%phi_P2_ref, this%gradt_P2_ref)
   call phi3_2d_mat (t_xy, this%phi_P3_ref, this%gradt_P3_ref)


   if (.not. is_curved ) then
      !  Rectilinear element
      call jacobian_p1_2d (t_xy, el_nds_xy, P2_NODES_PER_EL, xy_act, this%det, this%mat_B, this%mat_T, errco, emsg)
      call nberr%set(errco, emsg)
      RET_ON_NBERR_UNFOLD(nberr)

      if (this%det <= 0 .and. debug == 1) then
         write(*,*) "   !!!"
         write(*,*) "array_sol: det <= 0: i_el, det ", i_el, this%det
      endif

   else
      !  Isoparametric element, 2024-06-14 fix
      call jacobian_p2_2d (el_nds_xy, P2_NODES_PER_EL, this%phi_P2_ref, &
         this%gradt_P2_ref, xy_act, this%det, this%mat_B, this%mat_T, errco, emsg)
      call nberr%set(errco, emsg)
      RET_ON_NBERR_UNFOLD(nberr)

   endif


   !  grad_i  = gradient on the actual triangle
   !  grad_i  = Transpose(mat_T)*grad_i0
   !  Calculation of the matrix-matrix product:
   call DGEMM('Transpose','N', 2, 3,  2, D_ONE, this%mat_T, 2, &
      this%gradt_P1_ref, 2, D_ZERO, this%gradt_P1_act, 2)
   call DGEMM('Transpose','N', 2, 6,  2, D_ONE, this%mat_T, 2, &
      this%gradt_P2_ref, 2, D_ZERO, this%gradt_P2_act, 2)
   call DGEMM('Transpose','N', 2, 10, 2, D_ONE, this%mat_T, 2, &
      this%gradt_P3_ref, 2, D_ZERO, this%gradt_P3_act, 2)


end subroutine

! Evaluates the vector element and its transverse curl
! of basis function bf_j (1..3) of transverse entity ety_j (1..4)
! at the current position previously set by evaluate_at_position
subroutine  BasisFunctions_evaluate_vector_elts(this, bf_j, ety_j, vec_phi, curlt_phi)
   use numbatmod

   class(BasisFunctions) this
   integer(8) bf_j, ety_j
   double precision vec_phi(2)  ! gradient of a P
   double precision curlt_phi    ! \zhat \dot (\nabla_t \cross E_t)

   call evaluate_vector_elts(bf_j, ety_j, this%vector_elt_map, this%phi_P2_ref, &
      this%gradt_P1_act, this%gradt_P2_act, vec_phi, curlt_phi)

end subroutine
