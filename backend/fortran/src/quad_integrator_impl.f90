subroutine QuadIntegrator_setup_reference_quadratures(this)

   class(QuadIntegrator) this

   call quad_triangle (this%n_quad, NQUAD_MAX, &
      this%wt_quad, this%x_quad, this%y_quad)

end subroutine

subroutine QuadIntegrator_build_transforms_at(this, qi, nds_xy, &
   is_curved, errco, emsg)

   class(QuadIntegrator) this

   integer(8) qi  ! quad_index
   double precision nds_xy(2,P2_NODES_PER_EL)  ! positions of nodes of the element
   logical is_curved

   double precision xy_ref(2), xy_act(2), ww, det
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg


   ! locals
   xy_ref(1) = this%x_quad(qi)
   xy_ref(2) = this%y_quad(qi)
   ww = this%wt_quad(qi)

   !  xy_ref   = coordinate on the reference triangle
   !  xy_act = coordinate on the actual triangle

   !  phi_P2_ref = values of Lagrange polynomials (1-6) at each local node.
   !  gradt_P2_ref = gradient on the reference triangle (P2 element)
   call phi2_2d_mat(xy_ref, this%phi_P2_ref, this%gradt_P2_ref)

   if (.not. is_curved) then ! Rectilinear element
      call jacobian_p1_2d(xy_ref, nds_xy, P2_NODES_PER_EL, xy_act, &
         this%det, this%mat_B, this%mat_T, errco, emsg)
   else ! Isoparametric element !  fixed 2024/6/12
      call jacobian_p2_2d(nds_xy, P2_NODES_PER_EL, this%phi_P2_ref, this%gradt_P2_ref, xy_act, &
         this%det, this%mat_B, this%mat_T, errco, emsg)
   endif


   if(abs(this%det) .lt. 1.0d-20) then
      write(emsg,*) 'quadrature integration: det too small: ',  det
      errco = NBERR_BAD_QUAD_INT
      return
   endif


   ! gradt_P2_act  = gradient on the actual triangle
   ! gradt_P2_act  = Transpose(mat_T)*gradt_P2_ref
   ! Calculation of the matrix-matrix product:
   call DGEMM('Transpose','N', 2, 6, 2, D_ONE, this%mat_T, 2, &
      this%gradt_P2_ref, 2, D_ZERO, this%gradt_P2_act, 2)


end subroutine
