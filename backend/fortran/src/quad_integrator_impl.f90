subroutine QuadIntegrator_setup_reference_quadratures(this)

   class(QuadIntegrator) this

   call quad_triangle (this%n_quad, NQUAD_MAX, &
      this%wt_quad, this%x_quad, this%y_quad)

end subroutine

subroutine QuadIntegrator_build_transforms_at(this, qi, nds_xy, is_curved, do_P3, nberr)

   class(QuadIntegrator) this

   integer(8) qi  ! quad_index
   double precision nds_xy(2,P2_NODES_PER_EL)  ! positions of nodes of the element
   logical is_curved, do_P3

   double precision xy_ref(2), xy_act(2)
   type(NBError) nberr

   integer(8)  :: errco
   character(len=EMSG_LENGTH) :: emsg

   errco = 0

   xy_ref(1) = this%x_quad(qi)
   xy_ref(2) = this%y_quad(qi)

   !  xy_ref   = coordinate on the reference triangle
   !  xy_act = coordinate on the actual triangle

   !  phi_P2_ref = values of the P2 polynomials (1-6) at a given point in the reference triangle.
   !  gradt_P2_ref = gradients of the P2s at that point
   call phi2_2d_mat(xy_ref, this%phi_P2_ref, this%gradt_P2_ref)

   if (do_P3) then
      call phi3_2d_mat(xy_ref, this%phi_P3_ref, this%gradt_P3_ref)
   endif

   if (.not. is_curved) then ! Rectilinear element
      call jacobian_p1_2d(xy_ref, nds_xy, P2_NODES_PER_EL, xy_act, &
         this%det, this%mat_B, this%mat_T, errco, emsg)
   else ! Isoparametric element !  fixed 2024/6/12
      call jacobian_p2_2d(nds_xy, P2_NODES_PER_EL, this%phi_P2_ref, this%gradt_P2_ref, xy_act, &
         this%det, this%mat_B, this%mat_T, errco, emsg)
   endif

   call nberr%set(errco, emsg)
   RET_ON_NBERR_UNFOLD(nberr)

   if(abs(this%det) .lt. 1.0d-20) then
      write(emsg,*) 'quadrature integration: det too small: ',  this%det
      call nberr%set(NBERR_BAD_QUAD_INT, emsg)
      return
   endif

   ! gradt_P2_act  = gradient on the actual triangle
   ! gradt_P2_act  = Transpose(mat_T)*gradt_P2_ref
   ! Calculation of the matrix-matrix product:
   call DGEMM('Transpose','N', 2, 6, 2, D_ONE, this%mat_T, 2, &
      this%gradt_P2_ref, 2, D_ZERO, this%gradt_P2_act, 2)

   if (do_P3) then
      call DGEMM('Transpose','N', 2, 10, 2, D_ONE, this%mat_T, 2,this%gradt_P3_ref, 2, D_ZERO, this%gradt_P3_act, 2)
   endif

   this%quadweight = this%wt_quad(qi) * abs(this%det)

end subroutine

function QuadIntegrator_get_current_quadweight(this) result(res)
   class(QuadIntegrator) this

   double precision res        ! transformation determinant

   res = this%quadweight
end function


