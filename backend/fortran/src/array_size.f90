
#include "numbat_decl.h"

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Estimates the work space sizes that will be needed
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine array_size (n_msh_pts, n_msh_el, n_modes, &
 int_size, cmplx_size, real_size, n_ddl, errco, emsg)

      use numbatmod
      integer*8 n_msh_pts, n_msh_el
      integer*8 :: n_modes
      integer*8 int_size, cmplx_size, real_size, n_ddl

      integer errco
      character(len=EMSG_LENGTH) emsg


!     Local variables
      integer*8 nnodes, npt, npt_p3
      integer*8 nvect, ordre_ls
      integer*8 nonz, nonz_max, max_row_len
      integer*8 neq, neq_PW
      integer*8 n_edge, n_ddl_max
      integer*8 ltrav
!
!  Declare the pointers of the integer super-vector
      integer*8 ip_type_nod, ip_type_el, ip_table_nod
      integer*8 ip_table_E, ip_table_N_E_F, ip_visite
      integer*8 ip_type_N_E_F, ip_eq
      integer*8 ip_period_N, ip_nperiod_N
      integer*8 ip_period_N_E_F, ip_nperiod_N_E_F
      integer*8 ip_index_pw_inv
!  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_x_N_E_F, jp_rhs
!      integer*8 jp_matD, jp_matL, jp_matU
!      integer*8 jp_matD2, jp_matL2, jp_matU2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_eigen_modes_tmp, jp_trav, jp_vp, jp_eigen_pol
      integer*8 jp_overlap_L, jp_overlap_J, jp_overlap_J_dagger
      integer*8 jp_flux
      integer*8 jp_overlap_K, jp_X_mat
      integer*8 jp_sol1, jp_sol2, jp_sol1b
      integer*8 jp_sol1_H, jp_sol1b_H
      integer*8 jp_eigen_modes, jp_eigen_modes1, jp_eigen_modes2
      integer*8 jp_T, jp_R, jp_T12, jp_R12, jp_T21, jp_R21
      integer*8 jp_T_Lambda, jp_R_Lambda
      integer*8 jp_X_mat_b
!     Declare the pointers of the real super-vector
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im
      integer*8 kp_mat1_re, kp_mat1_im

!     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 jp_mat2
      integer*8 ip_work, ip_work_sort, ip_work_sort2


!f2py intent(in)  n_msh_el, n_modes

!f2py intent(out)  int_size, cmplx_size, real_size

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      nnodes = 6
      ordre_ls = 0
!
!     For most of the FEM meshes I have used, I have observed that:
!     npt is approximately equal to n_msh_el * 2.1
!
      npt = n_msh_el * 3
      n_edge = (npt + n_msh_el) / 2
      npt_p3 = npt + n_edge + n_msh_el
      nvect = 2*n_modes + n_modes/2 +3

!     Euler's polyhedron formula (no holes):
!     V - E + F = 2
!     V, E, and F are respectively the numbers of vertices (corners), edges and faces (triangles)
!
!     Since V - E = 2 - n_msh_el and V + E = npt, we have E = n_edge = (npt+n_msh_el-2)/2
!

      n_ddl = n_edge + n_msh_el + npt_p3
      neq = 3 * (n_edge + n_msh_el) + npt_p3
      neq_PW = (2*ordre_ls+1)**2
!     For most of the FEM meshes I have used, I have observed that:
!     nonz = 34.25 * neq
      nonz = 40 * neq
      nonz_max = nonz

!     I have observed that: max_row_len < 200
      max_row_len = 200
!ccccc
      ip_type_nod = 1
      ip_type_el = ip_type_nod + npt
!       ! pointer to FEM connectivity table
      ip_table_nod = ip_type_el + n_msh_el
      ip_table_N_E_F = ip_table_nod + nnodes*n_msh_el

      n_ddl_max = npt + n_msh_el
      ip_visite =  ip_table_N_E_F  + 14*n_msh_el
      ip_table_E = ip_visite + n_ddl_max

      ip_type_N_E_F = ip_table_E + 4*n_edge


        ip_period_N = ip_type_N_E_F + 2*n_ddl
        ip_nperiod_N = ip_period_N + npt
        ip_period_N_E_F = ip_nperiod_N + npt
        ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
        ip_eq = ip_nperiod_N_E_F + n_ddl

      ip_index_pw_inv = ip_eq + 3*n_ddl
      ip_col_ptr = ip_index_pw_inv + neq_PW
      ip_row = ip_col_ptr + neq + 1

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*n_ddl
      ip_work_sort2 = ip_work_sort + max_row_len

      int_size = ip_work_sort2 + max_row_len

!ccccc
      jp_x = 1
      jp_x_N_E_F = jp_x + 2*npt

      jp_rhs = jp_x_N_E_F + 3*n_ddl
!     jp_rhs will also be used (in gmsh_post_process) to store a solution
      jp_mat2 = jp_rhs + max(neq, 3*npt)

      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq

      jp_sol1 = jp_resid + neq
      jp_sol1b = jp_sol1 + 3*(nnodes+7)*n_modes*n_msh_el
      jp_sol2 = jp_sol1b + 3*(nnodes+7)*n_modes*n_msh_el
      jp_sol1_H = jp_sol2 + 3*(nnodes+7)*n_modes*n_msh_el
      jp_sol1b_H = jp_sol1_H + 3*nnodes*n_modes*n_msh_el
      jp_eigen_modes1 = jp_sol1b_H + 3*nnodes*n_modes*n_msh_el
      jp_eigen_modes2 = jp_eigen_modes1 + n_modes + 1
!       ! Eigenvectors
      jp_vschur = jp_eigen_modes2 + n_modes + 1
      jp_eigen_modes = jp_vschur + neq*nvect
      jp_eigen_pol = jp_eigen_modes + n_modes + 1
      jp_eigen_modes_tmp = jp_eigen_pol + n_modes*4
      jp_trav = jp_eigen_modes_tmp + n_modes + 1

      ltrav = 3*nvect*(nvect+2)
      jp_vp = jp_trav + ltrav
      jp_overlap_L = jp_vp + neq*n_modes
      jp_flux = jp_overlap_L + n_modes*n_modes
      jp_overlap_J = jp_flux + n_modes
      jp_overlap_J_dagger = jp_overlap_J + 2*neq_PW*n_modes
      jp_overlap_K = jp_overlap_J_dagger + 2*neq_PW*n_modes
      jp_X_mat = jp_overlap_K + 2*neq_PW*n_modes
      jp_T12 = jp_X_mat + 2*neq_PW*2*neq_PW
      jp_T21 = jp_T12 + 2*neq_PW*n_modes
      jp_R12 = jp_T21 + 2*neq_PW*n_modes
      jp_R21 = jp_R12 + 4*neq_PW*neq_PW
      jp_T = jp_R21 + n_modes*n_modes
      jp_R = jp_T + 2*neq_PW*2*neq_PW
      jp_T_Lambda = jp_R + 2*neq_PW*2*neq_PW
      jp_R_Lambda = jp_T_Lambda + 2*neq_PW*2*neq_PW

!c      if (substrate .eq. 1) then
      jp_X_mat_b = jp_R_Lambda + 2*neq_PW*2*neq_PW
      cmplx_size = jp_X_mat_b + 2*neq_PW*2*neq_PW
!c      else
!c        cmplx_size = jp_R_Lambda + 2*neq_PW*2*neq_PW
!c      endif

!ccccc

      kp_rhs_re = 1
      kp_rhs_im = kp_rhs_re + neq
      kp_lhs_re = kp_rhs_im + neq
      kp_lhs_im = kp_lhs_re + neq
      kp_mat1_re = kp_lhs_im + neq
      kp_mat1_im = kp_mat1_re + nonz
      real_size = kp_mat1_im + nonz

! cccccc
! c     SOME 32 bit integers for UMFPACK AND ARPACK
!       ip_col_ptr_32 = 1
!       ip_row_32 = ip_col_ptr_32 + neq + 1
!       int_size_32 = ip_row_32 + nonz

! cccccc
!       write(*,*) "array_size:"
!       write(*,*) "int_size = ", int_size
!       write(*,*) "cmplx_size = ", cmplx_size
!       write(*,*) "real_size = ", real_size
!       write(*,*) "int_size_32 = ", int_size_32

!       open (unit=26, file="Output/array_size.txt", status='unknown')
!         write(26,*) "int_size = ", int_size
!         write(26,*) "cmplx_size = ", cmplx_size
!         write(26,*) "real_size = ", real_size
!         write(26,*) "int_size_32 = ", int_size_32
!         write(26,*)
!         write(26,*) "n_msh_el = ", n_msh_el
!         write(26,*) "n_modes = ", n_modes
!         write(26,*) "nvect = ", nvect
!         write(26,*) "ordre_ls = ", ordre_ls
!         write(26,*)
!         write(26,*) "npt = ", npt
!         write(26,*) "n_edge = ", n_edge
!         write(26,*) "npt_p3 = ", npt_p3
!         write(26,*) "n_ddl = ", n_ddl
!         write(26,*) "neq = ", neq
!         write(26,*) "neq_PW = ", neq_PW
!         write(26,*) "nonz_max = ", nonz_max
!         write(26,*) "nonz = ", nonz
!         write(26,*) "max_row_len = ", max_row_len
!       close(26)
!
!



!     Dimension checking, can this actually happen?

      if ((3*n_msh_pts+n_msh_el+nnodes*n_msh_el) .gt. int_size) then
         write(emsg,*) "py_calc_modes(_AC): ", &
   "(3*n_msh_pts+n_msh_el+nnodes*n_msh_el) + ", &
   "n_msh_pts > int_max : ", &
    (3*n_msh_pts+n_msh_el+nnodes*n_msh_el), int_size, &
   "increase the size of int_max"
         errco = -2
         return
      endif

      if ((7*n_msh_pts) .gt. cmplx_size) then
         write(emsg,*) "py_calc_modes(_AC): (7*n_msh_pts)>cmplx_max : ", &
    (7*n_msh_pts), cmplx_size, &
   "increase the size of cmplx_max"
         errco = -2
         return
      endif


      end subroutine array_size

