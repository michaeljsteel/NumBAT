#include "numbat_decl.h"

subroutine AnalyticIntegrator_build_transforms_at(this, nds_xy, errco, emsg)

   class(AnalyticIntegrator) this

   double precision nds_xy(2,P2_NODES_PER_EL)
   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !  The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
   !  maps the current triangle to the reference triangle.

   this%mat_B(:, 1) = nds_xy(:, 2) - nds_xy(:, 1)
   this%mat_B(:, 2) = nds_xy(:, 3) - nds_xy(:, 1)

   this%det = this%mat_B(1,1) * this%mat_B(2,2) - this%mat_B(1,2) * this%mat_B(2,1)

   if (abs(this%det) .le. 1.0d-22) then
      write(emsg,*) 'analytic integrator: Determinant = 0 :', this%det, "nds_xy = ", nds_xy
      errco = NBERR_BAD_DETERMINANT
      return
   endif


   !  mat_T = mat_B ^{-1}
   this%mat_T(1,1) =   this%mat_B(2,2) / this%det
   this%mat_T(2,2) =   this%mat_B(1,1) / this%det
   this%mat_T(1,2) = - this%mat_B(1,2) / this%det
   this%mat_T(2,1) = - this%mat_B(2,1) / this%det

   !  Note that if grad_i_0 is the gradient on the reference triangle,
   !  then the gradient on the actual triangle is:
   !  grad_i  = Transpose(mat_T)*grad_i0
   !
   !  mat_T_tr = Transpose(mat_T) = Transpose(mat_B ^{-1})
   this%mat_T_tr(1,1) = this%mat_T(1,1)
   this%mat_T_tr(2,2) = this%mat_T(2,2)
   this%mat_T_tr(1,2) = this%mat_T(2,1)
   this%mat_T_tr(2,1) = this%mat_T(1,2)


end subroutine

subroutine PyFrontEnd_init_from_py(this, n_msh_el, n_msh_pts, elnd_to_mesh, v_nd_xy, errco, emsg)
   class(PyFrontEnd) this

   integer(8) n_msh_el
   integer(8) n_msh_pts
   integer(8) :: elnd_to_mesh(P2_NODES_PER_EL, n_msh_el)
   double precision:: v_nd_xy(2, n_msh_pts)

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   call integer_alloc_2d(this%elnd_to_mesh, P2_NODES_PER_EL, n_msh_el, 'elnd_to_mesh', errco, emsg);
   RETONERROR(errco)

   call double_alloc_2d(this%v_nd_xy, 2_8, n_msh_pts, 'v_nd_xy', errco, emsg);
   RETONERROR(errco)

   this%n_msh_pts = n_msh_pts
   this%n_msh_el = n_msh_el
   this%elnd_to_mesh = elnd_to_mesh
   this%v_nd_xy = v_nd_xy

end subroutine

subroutine PyFrontEnd_nodes_at_el(this, i_el, nds_xy)
   class(PyFrontEnd) this
   integer(8) i_el
   double precision nds_xy(2,P2_NODES_PER_EL)
   integer(8) j

   do j=1,P2_NODES_PER_EL
      nds_xy(:, j) = this%v_nd_xy(:, this%elnd_to_mesh(j,i_el))
   enddo

end subroutine
