#include "numbat_decl.h"
!  P2 basis function over the unit triangle

!  Quadratic basis function = P2 * Grad P1

!  vector_elt_map(1,j,i) = k : number on data to be stored: if k=3 only one gradient will be used; k=4 => 2 gradients
!  vector_elt_map(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  vector_elt_map(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial phi_n
!  vector_elt_map(4,j,i)     : it will be used only if k=4

!  The basis functions are stored in the 3rd index i as : centre, edges 1(4), 2(5), 3(6)
!                                 in the 2nd index j as 3 functions on each edge (for i=1), and along each edge (for i=2..4)


subroutine BasisFunctions_build_vector_elt_map(this, el_nds)
   use numbatmod

   class(BasisFunctions) this

   integer(8) el_nds(P2_NODES_PER_EL)     ! global indices of the P2 nodes at current element


   integer ety, bf, j
   integer(8) list_end(2,3), j2
   integer(8) ls_n(3)
   integer(8) edge_end, elo, ehi, iedge


!  Endpoints of the 6 edges (mid-point) of the reference triangle

   iedge = 1
   list_end(1,iedge) = 1
   list_end(2,iedge) = 2

   iedge = 2
   list_end(1,iedge) = 2
   list_end(2,iedge) = 3

   iedge = 3
   list_end(1,iedge) = 1
   list_end(2,iedge) = 3


   !  scan the element face
   ety=1

   !  The mid-edge nodes of the face
   !  scan the mid-edge nodes of the face

   do bf=1,3

      this%vector_elt_map(1,bf,ety) = 3        !  number on data to be stored

      this%vector_elt_map(2,bf,ety) = bf+3      !  the mid-edge number

      !  give the node opposite to the mid-edge node (j+3)
      j2 = modulo(bf+2,3)
      if( j2 .eq. 0 ) j2 = 3

      this%vector_elt_map(3,bf,ety) = j2
      ! second gradient not needed for this function
      this%vector_elt_map(4,bf,ety) = 0
   enddo

   !  scan the 3 element edges
   do ety=2,4
      !  2 end-point basis vectors are attached to the edge i
      !  scan the end nodes of the edge

      iedge = ety-1

      ! Find the indices corresponding to the nodes at the ends of the current edge
      ! They need to be sorted by the absolute value of the nodes in the overal mesh table (via el_nds)
      do j=1,2
         edge_end = list_end(j,iedge)
         ls_n(j) = el_nds(edge_end)
      enddo

      elo=1
      ehi=2
      ! swap if absolute node ordering is reversed
      if (ls_n(1) .gt. ls_n(2)) then
         elo=2
         ehi=1
      endif

      elo = list_end(elo,iedge)
      ehi = list_end(ehi,iedge)

      bf = 1
      this%vector_elt_map(1,bf,ety) = 3
      this%vector_elt_map(2,bf,ety) = elo
      this%vector_elt_map(3,bf,ety) = ehi
      this%vector_elt_map(4,bf,ety) = 0

      bf = 2
      this%vector_elt_map(1,bf,ety) = 3
      this%vector_elt_map(2,bf,ety) = ehi
      this%vector_elt_map(3,bf,ety) = elo
      this%vector_elt_map(4,bf,ety) = 0

      bf = 3
      this%vector_elt_map(1,bf,ety) = 4
      this%vector_elt_map(2,bf,ety) = iedge+3   !  add 3 to get the correct edge number
      this%vector_elt_map(3,bf,ety) = elo
      this%vector_elt_map(4,bf,ety) = ehi

   enddo

end subroutine

subroutine BasisFunctions_set_affine_for_elt(this, el_nds_xy, nberr)
   use numbatmod

   class(BasisFunctions) this
   double precision el_nds_xy(2,P2_NODES_PER_EL)
   type(NBError) nberr

   !double precision t_xy(2), xy_act
   character(len=EMSG_LENGTH) :: emsg


   this%mat_B(:,1) = el_nds_xy(:,2) - el_nds_xy(:,1)
   this%mat_B(:,2) = el_nds_xy(:,3) - el_nds_xy(:,1)

   this%det = this%mat_B(1,1) * this%mat_B(2,2) - this%mat_B(1,2) * this%mat_B(2,1)


   if (abs(this%det) .le. 1.0d-22) then
      write(emsg,*) 'Bad jacobian in set_affine_for_elt'
      call nberr%set(NBERR_BAD_JACOBIAN, emsg)
      return

   endif

   !   mat_T = Inverse of mat_B
   this%mat_T(1,1) =  this%mat_B(2,2) / this%det
   this%mat_T(2,2) =  this%mat_B(1,1) / this%det
   this%mat_T(1,2) = -this%mat_B(1,2) / this%det
   this%mat_T(2,1) = -this%mat_B(2,1) / this%det

   ! !	mat_T_tr = Tanspose(mat_T)
   this%mat_T_tr(1,1) = this%mat_T(1,1)
   this%mat_T_tr(1,2) = this%mat_T(2,1)
   this%mat_T_tr(2,1) = this%mat_T(1,2)
   this%mat_T_tr(2,2) = this%mat_T(2,2)

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
      RET_ON_NBERR(nberr)

      if (this%det <= 0 .and. debug == 1) then
         write(*,*) "   !!!"
         write(*,*) "array_sol: det <= 0: i_el, det ", i_el, this%det
      endif

   else
      !  Isoparametric element, 2024-06-14 fix
      call jacobian_p2_2d (el_nds_xy, P2_NODES_PER_EL, this%phi_P2_ref, &
         this%gradt_P2_ref, xy_act, this%det, this%mat_B, this%mat_T, errco, emsg)
      call nberr%set(errco, emsg)
      RET_ON_NBERR(nberr)

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

!  P2 vector basis functions over the reference unit triangle

!  Compute:
!  a quadratic vector basis function (vec_phi = P2 * Grad P1) and its transverse curl (curlt_phi)

! vector_elt_map indexes the pieces of the basis function
!  vector_elt_map(1,j,i) = k : number of gradients needed: if k=3 only one gradient will be used; k=4 => 2 gradients
!  vector_elt_map(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
!  vector_elt_map(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial eta_n
!  vector_elt_map(4,j,i)     : gradient of the second P1 polynomial eta_n2 if needed


! Generate the vector basis function for the bf_j bf_j of transverse entity ety_trans and its transverse curl (purely z-component)
! in terms of the P1 and P2 scalar functions

! formulas
! ety_trans=1, face element
!    3 edge functions at P2 nodes 4, 5, 6
!       j=1, node 4: phi_4 Grad eta_3
!       j=2, node 5: phi_5 Grad eta_1
!       j=3, node 6: phi_6 Grad eta_2

! ety_trans=2, edge 1, end nodes 1&2
!       j=1, : phi_1 Grad phi_2
!       j=2, : phi_2 Grad phi_1
!       j=3, : phi_4 (Grad phi_1-Grad phi_2)

! ety_trans=3, edge 2, end nodes 2&3
!       j=1, : phi_2 Grad phi_3
!       j=2, : phi_3 Grad phi_2
!       j=3, : phi_5 (Grad phi_2-Grad phi_3)

! ety_trans=4, edge 3, end nodes 3&1
!       j=1, : phi_3 Grad phi_1
!       j=2, : phi_1 Grad phi_3
!       j=3, : phi_6 (Grad phi_1-Grad phi_3)

subroutine  BasisFunctions_evaluate_vector_elts(this, bf_j, ety_trans, vec_phi, curlt_phi)
   use numbatmod

   class(BasisFunctions) this
   integer(8) bf_j, ety_trans
   double precision vec_phi(2)  ! gradient of a P
   double precision curlt_phi    ! \zhat \dot (\nabla_t \cross E_t)

   double precision grad_p1(2), grad_p2(2), phi
   integer(8) k, m, n1, n2

   ! call evaluate_vector_elts(bf_j, ety_j, this%vector_elt_map, this%phi_P2_ref, &
   !    this%gradt_P1_act, this%gradt_P2_act, vec_phi, curlt_phi)


   k  = this%vector_elt_map(1, bf_j, ety_trans)
   m  = this%vector_elt_map(2, bf_j, ety_trans)
   n1 = this%vector_elt_map(3, bf_j, ety_trans)
   n2 = this%vector_elt_map(4, bf_j, ety_trans)



   if (k .eq. 3) then

      phi = this%phi_P2_ref(m)
      grad_p2 = this%gradt_P2_act(:,m)
      grad_p1 = this%gradt_P1_act(:,n1)
      vec_phi = phi * grad_p1


      !elseif (k .eq. 4) then
   else ! k==4

      phi = this%phi_P2_ref(m)
      grad_p2 = this%gradt_P2_act(:,m)
      grad_p1 = this%gradt_P1_act(:,n1) - this%gradt_P1_act(:,n2)
      vec_phi = phi * grad_p1


   endif

   !  Curl_t E = Det( grad_p2,  grad_p1)
   curlt_phi = grad_p2(1)*grad_p1(2) - grad_p2(2)*grad_p1(1)


end subroutine

subroutine BasisFunctions_find_derivatives(this, idof, ety, &
   vec_phi_x, curlt_phi_x, phi_P3_x, gradt_P3_x)

   use numbatmod
   class(BasisFunctions) this

   integer(8) idof, ety

   double precision vec_phi_x(2), curlt_phi_x
   double precision phi_P3_x
   double precision gradt_P3_x(2)

   if (ety .le. N_ETY_TRANSVERSE) then ! A transverse dof (edge or face)
      ! Uses P2 vector elements so determine the basis vector
      call this%evaluate_vector_elts(idof, ety, vec_phi_x, curlt_phi_x)

      phi_P3_x = D_ZERO
      gradt_P3_x = D_ZERO
   else   ! a longitudinal dof, use P3 scalar element
      vec_phi_x = D_ZERO
      curlt_phi_x = D_ZERO

      phi_P3_x = this%phi_P3_ref(ety-N_ETY_TRANSVERSE)
      gradt_P3_x(:) = this%gradt_P3_act(:,ety-N_ETY_TRANSVERSE)
   endif

end subroutine



subroutine BasisFunctions_get_triint_p2_p2(this, m_p2_p2)

   class(BasisFunctions) this
   double precision m_p2_p2(P2_NODES_PER_EL, P2_NODES_PER_EL)

   call find_overlaps_p2_p2(m_p2_p2, this%det)

end subroutine

subroutine BasisFunctions_get_triint_set_p2_p2(this)
   class(BasisFunctions) this
    call find_overlaps_p2_p2(this%p2_p2, this%det)
    call find_overlaps_p2_p2x(this%p2_p2x, this%mat_T_tr, this%det)
    call find_overlaps_p2_p2y(this%p2_p2y, this%mat_T_tr, this%det)
    call find_overlaps_p2x_p2x(this%p2x_p2x, this%mat_T_tr, this%det)
    call find_overlaps_p2x_p2y(this%p2x_p2y, this%mat_T_tr, this%det)
    call find_overlaps_p2y_p2y(this%p2y_p2y, this%mat_T_tr, this%det)
end subroutine
