#include "numbat_decl.h"


subroutine calc_AC_modes(n_modes, q_ac, dimscale_in_m, shift_nu, &
   i_bnd_cdns, itermax, tol, debug, show_mem_est, &
   symmetry_flag, n_typ_el, &
   c_tensor, rho, &
   supplied_geo_flag, &
   mesh_file, n_msh_pts, n_msh_el, &
   type_nod, &
   table_nod, type_el, mesh_xy, &
   v_eigs_nu, sol1, mode_pol, errco, emsg)

   ! q_ac :   acoustic wave number (q_ac)
   ! n_modes:  desired number of solved acoustic modes
   ! n_msh_pts:  number of nodes in mesh
   ! n_msh_el:   number of (triang) elements in mesh
   ! n_type_el:  number of types of material
   ! type_nod:    ??
   ! table_nod:
   ! type_el:
   ! mesh_xy
   ! v_eigs_nu:  eigen frequencies nu=omega/(2D_PI) for each mode
   ! sol1:
   ! mode_pol:

   !***********************************************************************
   !
   ! Program:
   ! FEM solver of Acoustic waveguide problems.
   ! This subroutine is compiled by f2py & called in mode_calcs.py
   !
   ! Authors:
   ! Bjorn Sturmberg & Kokou B. Dossou
   !
   !***********************************************************************
   !
   use numbatmod
   !f2py integer, parameter :: EMSG_LENGTH
   integer*8 :: n_modes

   complex*16 q_ac
   double precision dimscale_in_m
   complex*16 shift_nu
   integer*8 int_max, cmplx_max, int_used, cmplx_used
   integer*8 real_max, real_used, n_ddl
   integer errco
   character(len=EMSG_LENGTH) emsg



   integer*8 n_msh_el, n_msh_pts, i_bnd_cdns, neq
   integer*8 debug, show_mem_est

   !     ! Number of nodes per element
   integer*8 , parameter :: nodes_per_el =6
   integer*8 type_nod(n_msh_pts), type_el(n_msh_el)
   integer*8 table_nod(nodes_per_el, n_msh_el)

   double precision mesh_xy(2,n_msh_pts)
   complex*16, target :: sol1(3,nodes_per_el,n_modes,n_msh_el)
   complex*16, target :: v_eigs_nu(n_modes)
   complex*16 mode_pol(4,n_modes)

   integer*8 supplied_geo_flag, symmetry_flag



   integer :: alloc_stat

   integer*8, dimension(:), allocatable :: a_iwork
   complex*16, dimension(:), allocatable :: b_zwork
   double precision, dimension(:), allocatable :: c_dwork

   !double precision, dimension(:,:), allocatable :: d_dwork
   !complex*16, dimension(:,:), allocatable :: dummy_overlap_L  ! not actually used.

   !  Declare the pointers of the integer super-vector
   integer*8 ip_eq
   integer*8 ip_visite

   !  Declare the pointers of the real super-vector
   integer*8 jp_x, jp_mat2
   integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
   integer*8 jp_trav, jp_vp, jp_rhs
   integer*8 jp_eigenum_modes_tmp, jp_eigen_pol

   !     Declare the pointers of the real super-vector
   integer*8 kp_mat1_re, kp_mat1_im
   integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im

   !     Declare the pointers of for sparse matrix storage
   integer*8 ip_col_ptr, ip_row
   integer*8 ip_work, ip_work_sort, ip_work_sort2
   integer*8 nonz, nonz_max, max_row_len

   integer*8 n_typ_el
   complex*16 c_tensor(6,6,n_typ_el)
   complex*16 rho(n_typ_el)

   integer*8 i, j, ip
   integer*8 ui,  namelength


   double precision lat_vecs(2,2)
   double precision dim_x, dim_y

   complex*16 shift_omsq
   integer*8  i_base


   !  Variable used by valpr
   integer*8 ltrav, n_conv
   complex*16 z_beta, z_tmp, z_tmp0
   integer*8, dimension(:), allocatable :: iindex
   !     variable used by UMFPACK

   double precision time1, time2
   double precision stime1, stime2
   character*(8) start_date, end_date
   character*(10) start_time, end_time

   !  Variable used by valpr
   integer*8  nvect, itermax
   double precision tol

   !  Names and Controls

   character(len=FNAME_LENGTH)  mesh_file, gmsh_file, log_file, gmsh_file_pos


   integer :: is_em


   !f2py intent(in) q_ac, n_modes
   !f2py intent(in) debug, mesh_file, n_msh_pts, n_msh_el
   !f2py intent(in) dimscale_in_m, shift
   !f2py intent(in) i_bnd_cdns, itermax, tol
   !f2py intent(in) plot_modes, c_tensor, rho
   !f2py intent(in) cmplx_max, real_max, int_max
   !f2py intent(in) n_typ_el, supplied_geo_flag
   !f2py intent(in) type_nod, table_nod, type_el, mesh_xy, symmetry_flag

   !  Note: the dependent variables must be listed AFTER the
   !  independent variables that they depend on in the function call!

   !f2py depend(c_tensor) n_typ_el
   !f2py depend(rho) n_typ_el
   !f2py depend(type_nod) n_msh_pts
   !f2py depend(table_nod) nodes_per_el, n_msh_el
   !f2py depend(type_el) n_msh_el
   !f2py depend(mesh_xy) n_msh_pts

   !f2py intent(out) v_eigs_nu
   !f2py intent(out) sol1, mode_pol, table_nod, type_el, mesh_xy

   !f2py intent(out) errco
   !f2py intent(out) emsg

   !
   !CCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCC
   !
   !     Set parameter for the super-vectors of integer and real numbers
   !
   !       !ui = Unite dImpression
   ui = 6
   !      nodes_per_el = 6 ! Number of nodes per element


   !       nvect = 2*n_modes + n_modes/2 +3
   nvect = 3*n_modes + 3
   !
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   errco= 0
   emsg = ""
   debug=1


   call array_size(n_msh_pts, n_msh_el, n_modes, int_max, cmplx_max, real_max, &
      n_ddl, errco, emsg)
   RETONERROR(errco)

   allocate(a_iwork(int_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, int_max, "a", -1, errco, emsg)
   RETONERROR(errco)

   allocate(b_zwork(cmplx_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, cmplx_max, "b", -1, errco, emsg)
   RETONERROR(errco)

   allocate(c_dwork(real_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, real_max, "c", -1, errco, emsg)
   RETONERROR(errco)

   allocate(iindex(n_modes), STAT=alloc_stat)
   call check_alloc(alloc_stat, n_modes, "iindex", -1, errco, emsg)
   RETONERROR(errco)

   is_em = 0
   !    call prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, int_max, cmplx_max, real_max, &
   !       a_iwork, b_zwork, c_dwork, d_dwork, iindex, dummy_overlap_L, errco, emsg)
   !    RETONERROR(errco)

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !
   !     clean mesh_format
   namelength = len_trim(mesh_file)
   gmsh_file = mesh_file(1:namelength-5)//'.msh'
   gmsh_file_pos = mesh_file(1:namelength)
   log_file = mesh_file(1:namelength-5)//'-AC.log'
   if (debug .eq. 1) then
      write(*,*) "mesh_file = ", mesh_file
      write(*,*) "gmsh_file = ", gmsh_file
   endif

   !       ! initial time  in unit = sec.
   call cpu_time(time1)
   call date_and_time ( start_date, start_time )
   !
   !      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)

   dim_x = dimscale_in_m
   dim_y= dimscale_in_m
   shift_omsq= (2*D_PI*shift_nu)**2

   !####################  Start FEM PRE-PROCESSING  #######################

   !       ! pointer to FEM connectivity table
   ip_visite = 1
   ip_eq = ip_visite + n_msh_pts
   jp_x = 1
   !
   if (supplied_geo_flag .eq. 0) then
      call geometry (n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, dim_x, dim_y,  &
         mesh_file, mesh_xy, type_nod, type_el, table_nod, errco, emsg)
      if (errco .ne. 0) then
         return
      endif
   endif

   call lattice_vec (n_msh_pts, mesh_xy, lat_vecs, debug)

   !       if (debug .eq. 1) then
   !       open (unit=64, file="msh_check.txt",
   !      *         alloc_status="unknown")
   !         do i=1,n_msh_el
   !           write(64,*) i, type_el(i)
   !         enddo
   !         write(64,*)
   !         write(64,*)
   !         write(64,*)
   !         do i=1,n_msh_el
   !           do j=1,nodes_per_el
   !             write(64,*) i, j, table_nod(j,i)
   !           enddo
   !         enddo
   !         write(64,*)
   !         write(64,*)
   !         write(64,*)
   !         do j=1,nodes_per_el
   !           write(64,*) j, type_nod(j)
   !         enddo
   !       close(63)
   !       endif



   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes_AC: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
   endif

   !     Determine number of boundary conditions (neq) and 2D index array
   !     a_iwork(ip_eq)
   call bound_cond_AC (i_bnd_cdns, n_msh_pts, neq, type_nod, a_iwork(ip_eq), debug)
   !
   !
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !
   !     Sparse matrix storage
   ip_col_ptr = ip_eq + 3*n_msh_pts

   call csr_max_length_AC (n_msh_el, n_msh_pts, neq, nodes_per_el, &
      table_nod, a_iwork(ip_eq), a_iwork(ip_col_ptr), nonz_max)

   ip = ip_col_ptr + neq + 1
   if (ip .gt. int_max) then
      write(emsg,*) "py_calc_modes_AC: ip > int_max : ", &
         ip, int_max, &
         "py_calc_modes_AC: nonz_max = ", nonz_max, &
         "py_calc_modes_AC: increase the size of int_max"
      errco = -3
      return
   endif
   !
   ip_row = ip_col_ptr + neq + 1

   call csr_length_AC (n_msh_el, n_msh_pts, neq, nodes_per_el, &
      table_nod, a_iwork(ip_eq), a_iwork(ip_row), a_iwork(ip_col_ptr), nonz_max, &
      nonz, max_row_len, ip, int_max, debug)

   ip_work = ip_row + nonz
   ip_work_sort = ip_work + 3*n_msh_pts
   ip_work_sort2 = ip_work_sort + max_row_len

   !     sorting csr ...
   call sort_csr (neq, nonz, max_row_len, a_iwork(ip_row), &
      a_iwork(ip_col_ptr), a_iwork(ip_work_sort), a_iwork(ip_work), &
      a_iwork(ip_work_sort2))

   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes_AC: nonz_max = ", nonz_max
      write(ui,*) "py_calc_modes_AC: nonz = ", nonz
      write(ui,*) "py_calc_modes_AC: cmplx_max/nonz = ", &
         dble(cmplx_max)/dble(nonz)
   endif

   int_used = ip_work_sort2 + max_row_len

   if (int_max .lt. int_used) then
      write(emsg,*)"The size of the integer supervector is too small", &
         "integer super-vec: int_max  = ", int_max, &
         "integer super-vec: int_used = ", int_used
      errco = -4
      return
   endif

   jp_rhs = jp_x + 2*n_msh_pts
   !     jp_rhs will also be used (in gmsh_post_process) to store a solution
   jp_mat2 = jp_rhs + max(neq, 3*n_msh_pts)
   jp_vect1 = jp_mat2 + nonz
   jp_vect2 = jp_vect1 + neq
   jp_workd = jp_vect2 + neq
   jp_resid = jp_workd + 3*neq
   jp_eigenum_modes_tmp = jp_resid+3*nodes_per_el*n_modes*n_msh_el
   !       ! Eigenvectors
   jp_vschur = jp_eigenum_modes_tmp + n_modes + 1
   jp_eigen_pol = jp_vschur + neq*nvect
   jp_trav = jp_eigen_pol + n_modes*4

   ltrav = 3*nvect*(nvect+2)
   jp_vp = jp_trav + ltrav
   cmplx_used = jp_vp + neq*n_modes

   if (cmplx_max .lt. cmplx_used)  then
      write(emsg,*)"The size of the complex supervector is too small", &
         "complex super-vec: cmplx_max  = ", cmplx_max, &
         "complex super-vec: cmplx_used = ", cmplx_used
      errco = -5
      return
   endif

   !
   !ccccccccccccccccccccccccccccccccccccccccccccccccc
   !
   kp_rhs_re = 1
   kp_rhs_im = kp_rhs_re + neq
   kp_lhs_re = kp_rhs_im + neq
   kp_lhs_im = kp_lhs_re + neq
   kp_mat1_re = kp_lhs_im + neq
   kp_mat1_im = kp_mat1_re + nonz
   real_used = kp_mat1_im + nonz

   if (real_max .lt. real_used) then
      write(ui,*)
      write(ui,*) "The size of the real supervector is too small"
      write(ui,*) "2*nonz  = ", 2*nonz
      write(ui,*) "real super-vec: real_max  = ", real_max
      write(ui,*) "real super-vec: real_used = ", real_used

      write(emsg,*)"The size of the real supervector is too small", &
         "2*nonz  = ", 2*nonz, &
         "real super-vec: real_max  = ", real_max, &
         "real super-vec: real_used = ", real_used

      errco = -6
      return
   endif

   !
   !###############################################
   !
   !       ----------------------------------------------------------------
   !       convert from 1-based to 0-based
   !       ----------------------------------------------------------------
   !
   do 60 j = 1, neq+1
      a_iwork(j+ip_col_ptr-1) = a_iwork(j+ip_col_ptr-1) - 1
60 continue
   do 70 j = 1, nonz
      a_iwork(j+ip_row-1) = a_iwork(j+ip_row-1) - 1
70 continue
   !
   !
   !     The CSC iindexing, i.e., ip_col_ptr, is 1-based
   !       (but valpr.f will change the CSC iindexing to 0-based iindexing)
   i_base = 0

   !#####################  End FEM PRE-PROCESSING  #########################
   !
   write(ui,*)
   write(ui,*) "-----------------------------------------------"
   !       write(ui,*) " AC FEM, k_AC : ", real(q_ac), " 1/m"
   !       write(ui,*) "-----------------------------------------------"
   !       write(ui,*)

   !       if (debug .eq. 1) then
   !         write(ui,*) "py_calc_modes_AC: call to asmbly"
   !       endif
   write(ui,*) "AC FEM: "
   write(ui,*) "      - assembling linear system"
   !     Assemble the coefficient matrix K and M of the finite element equations

   call get_clocks(stime1, time1)

   call asmbly_AC (i_base, n_msh_el, n_msh_pts, neq, nodes_per_el, &
      shift_omsq, q_ac, n_typ_el, rho, c_tensor, &
      table_nod, type_el, a_iwork(ip_eq), &
      mesh_xy, nonz, a_iwork(ip_row), a_iwork(ip_col_ptr), &
      c_dwork(kp_mat1_re), c_dwork(kp_mat1_im), b_zwork(jp_mat2), a_iwork(ip_work), &
      symmetry_flag, debug)

   call get_clocks(stime2, time2)
   write(ui,'(A,F6.2,A)') '           cpu time  = ', (time2-time1), &
      ' secs.'
   write(ui,'(A,F6.2,A)') '           wall time = ', (stime2-stime1), &
      ' secs.'

   !       if (debug .eq. 1) then
   !         write(ui,*) "py_calc_modes_AC: call to valpr"
   !       endif
   write(ui,*) "      - solving linear system"
   call get_clocks(stime1, time1)


   call valpr_64_AC (i_base, nvect, n_modes, neq, itermax, ltrav, &
      tol, nonz, a_iwork(ip_row), a_iwork(ip_col_ptr), c_dwork(kp_mat1_re), &
      c_dwork(kp_mat1_im), b_zwork(jp_mat2), &
      b_zwork(jp_vect1), b_zwork(jp_vect2), b_zwork(jp_workd), b_zwork(jp_resid), &
      b_zwork(jp_vschur), v_eigs_nu, b_zwork(jp_trav), b_zwork(jp_vp), &
      c_dwork(kp_rhs_re), c_dwork(kp_rhs_im), c_dwork(kp_lhs_re), c_dwork(kp_lhs_im), n_conv, &
      debug, show_mem_est, errco, emsg)


   if (errco .ne. 0) then
      return
   endif

   call get_clocks(stime2, time2)
   write(ui,'(A,F6.2,A)') '           cpu time  = ', (time2-time1), &
      ' secs.'
   write(ui,'(A,F6.2,A)') '           wall time = ', (stime2-stime1), &
      ' secs.'


   if (n_conv .ne. n_modes) then
      write(emsg, '(A,I5,I5)') &
         "py_calc_modes_AC: convergence problem " // &
         "in valpr_64: n_conv != n_modes  ", n_conv, n_modes
      errco = -7
      return
   endif
   !
   do i=1,n_modes
      z_tmp0 = v_eigs_nu(i)
      z_tmp = 1.0d0/z_tmp0+shift_omsq
      z_beta = sqrt(z_tmp) / (2.0d0 * D_PI)
      !       Frequency (z_beta) should always be positive.
      if (dble(z_beta) .lt. 0) z_beta = -z_beta
      v_eigs_nu(i) = z_beta
   enddo
   !
   call z_indexx_AC (n_modes, v_eigs_nu, iindex)
   !
   !       The eigenvectors will be stored in the array sol1
   !       The eigenum_modesues and eigenvectors will be renumbered
   !                 using the permutation vector iindex
   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes_AC: call to array_sol"
   endif
   call array_sol_AC (n_modes, n_msh_el, n_msh_pts, neq, &
      nodes_per_el, iindex, table_nod, type_el, a_iwork(ip_eq), mesh_xy, &
      v_eigs_nu,  b_zwork(jp_eigenum_modes_tmp), mode_pol, b_zwork(jp_vp), sol1)

   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes_AC: array_sol returns call"
   endif
   !
   if(debug .eq. 1) then
      write(ui,*) 'iindex = ', (iindex(i), i=1,n_modes)
   endif
   if(debug .eq. 1) then
      write(ui,*)
      !         write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
      !         write(ui,*) "sqrt(shift_omsq)/(2*D_PI) = ", sqrt(omsq) / (2.0d0 * D_PI)
      do i=1,n_modes
         write(ui,"(i4,2(g22.14),2(g18.10))") i, v_eigs_nu(i)
      enddo
   endif

   !C    Save Original solution
   !      if (plot_modes .eq. 1) then
   !        dir_name = "AC_fields"
   !C        call write_sol_AC (n_modes, n_msh_el, nodes_per_el, lambda,
   !C      *       v_eigs_nu, sol1, mesh_file, dir_name)
   !C        call write_param (lambda, n_msh_pts, n_msh_el, i_bnd_cdns,
   !C    *       n_modes, nvect, itermax, tol, shift_omsq, lx, ly,
   !C    *       mesh_file, n_conv, dir_name)
   !        tchar = "AC_fields/All_plots_png_abs2_eE.geo"
   !        open (unit=34,file=tchar)
   !          do i=1,n_modes
   !            call gmsh_post_process_AC (i, n_modes, n_msh_el,
   !     *         n_msh_pts, nodes_per_el, table_nod, type_el,
   !     *         mesh_xy, v_eigs_nu, sol1, b_zwork(jp_rhs), a_iwork(ip_visite),
   !     *         gmsh_file_pos, dir_name, dimscale_in_m, debug)
   !          enddo
   !        close (unit=34)
   !      endif
   !C
   !
   !#########################  End Calculations  ###########################
   !
   call date_and_time ( end_date, end_time )
   call cpu_time(time2)
   !
   if (debug .eq. 1) then
      write(ui,*)
      write(ui,*) 'Total CPU time (sec.)  = ', (time2-time1)
      !
      open (unit=26,file=log_file)
      write(26,*)
      write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
      write(26,*) "Start date and time   = ", start_date, &
         " ; ", start_time
      write(26,*) "End date and time     = ", end_date, &
         " ; ", end_time
      write(26,*) "Total CPU time (sec.) = ",  (time2-time1)
      write(26,*)
      write(26,*) "q_ac = ", q_ac
      write(26,*) "shift_omsq= ", shift_omsq
      write(26,*)
      write(26,*) "n_msh_pts, n_msh_el, nodes_per_el  = ", n_msh_pts, &
         n_msh_el, nodes_per_el
      write(26,*) "neq, i_bnd_cdns = ", neq, i_bnd_cdns
      write(26,*) " lat_vecs:  = "
      write(26,"(2(f18.10))") lat_vecs
      write(26,*) "mesh_file = ", mesh_file
      write(26,*) "gmsh_file = ", gmsh_file
      write(26,*) "log_file  = ", log_file
      close(26)
      !
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      . (d=",dimscale_in_m,")"
      write(ui,*) "  and   we're   done!"
   endif

   write(ui,*) "-----------------------------------------------"
   write(ui,*)
   !
   deallocate(a_iwork, b_zwork, c_dwork, iindex)

end subroutine calc_AC_modes
