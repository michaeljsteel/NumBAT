#include "numbat_decl.h"

subroutine AnalyticIntegrator_build_transforms_at(this, nds_xy, nberr)

   class(AnalyticIntegrator) this

   double precision nds_xy(2,P2_NODES_PER_EL)
   type(NBError), intent(out) :: nberr

   character(len=EMSG_LENGTH) :: emsg

   ! The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
   ! maps the current triangle to the reference triangle.

   this%mat_B(:, 1) = nds_xy(:, 2) - nds_xy(:, 1)
   this%mat_B(:, 2) = nds_xy(:, 3) - nds_xy(:, 1)

   this%det = this%mat_B(1,1) * this%mat_B(2,2) - this%mat_B(1,2) * this%mat_B(2,1)

   if (abs(this%det) .le. 1.0d-20) then
      write(emsg,*) 'analytic integrator: Determinant = 0 :', this%det, "nds_xy = ", nds_xy
      call nberr%set(NBERR_BAD_DETERMINANT, emsg)
      return
   endif


   ! mat_T = mat_B ^{-1}
   this%mat_T(1,1) =   this%mat_B(2,2) / this%det
   this%mat_T(2,2) =   this%mat_B(1,1) / this%det
   this%mat_T(1,2) = - this%mat_B(1,2) / this%det
   this%mat_T(2,1) = - this%mat_B(2,1) / this%det

   ! Note that if grad_i_0 is the gradient on the reference triangle,
   ! then the gradient on the actual triangle is:
   ! grad_i  = Transpose(mat_T)*grad_i0
   !
   ! mat_T_tr = Transpose(mat_T) = Transpose(mat_B ^{-1})
   this%mat_T_tr(1,1) = this%mat_T(1,1)
   this%mat_T_tr(2,2) = this%mat_T(2,2)
   this%mat_T_tr(1,2) = this%mat_T(2,1)
   this%mat_T_tr(2,1) = this%mat_T(1,2)


end subroutine


subroutine PyFrontEnd_init_from_py(this, n_msh_elts, n_msh_pts, elnd_to_mshpt, v_mshpt_xy, nberr)

   class(PyFrontEnd) this
   type(NBError), intent(out) :: nberr

   integer(8) n_msh_elts
   integer(8) n_msh_pts
   integer(8) :: elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision:: v_mshpt_xy(2, n_msh_pts)


   call nberr%reset()

   call this%init_from_py2(n_msh_elts, n_msh_pts, elnd_to_mshpt, v_mshpt_xy, nberr)

end subroutine

subroutine PyFrontEnd_init_from_py2(this, n_msh_elts, n_msh_pts, elnd_to_mshpt, v_mshpt_xy, nberr)

   use numbatmod
   class(PyFrontEnd) this

   type(NBError) nberr
   integer(8) n_msh_elts
   integer(8) n_msh_pts
   integer(8) :: elnd_to_mshpt(P2_NODES_PER_EL, n_msh_elts)
   double precision:: v_mshpt_xy(2, n_msh_pts)

   call integer_alloc_2d(this%elnd_to_mshpt, P2_NODES_PER_EL, n_msh_elts, 'elnd_to_mshpt', nberr);
   RET_ON_NBERR(nberr)

   call double_alloc_2d(this%v_mshpt_xy, 2_8, n_msh_pts, 'v_mshpt_xy', nberr);
   RET_ON_NBERR(nberr)

   this%n_msh_pts = n_msh_pts
   this%n_msh_elts = n_msh_elts
   this%elnd_to_mshpt = elnd_to_mshpt
   this%v_mshpt_xy = v_mshpt_xy

end subroutine

subroutine PyFrontEnd_nodes_at_el(this, i_el, nds_xy)
   class(PyFrontEnd) this
   integer(8) i_el
   double precision nds_xy(2,P2_NODES_PER_EL)
   integer(8) j


   do j=1,P2_NODES_PER_EL
      nds_xy(:, j) = this%v_mshpt_xy(:, this%elnd_to_mshpt(j,i_el))
   enddo

   this%cur_elt = i_el
   this%cur_nds_xy = nds_xy

end subroutine

function PyFrontEnd_elt_is_curved(this) result(iscurv)
   class(PyFrontEnd) this
   logical iscurv

   iscurv = log_is_curved_elem_tri (P2_NODES_PER_EL, this%cur_nds_xy)
end function
