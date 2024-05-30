#include "numbat_decl.h"


      subroutine calc_AC_modes(
c     Explicit inputs
     *    num_modes, q_ac, 
     *    d_in_m, shift_nu,
     *    mesh_file, n_msh_pts, n_msh_el,
     *    symmetry_flag, n_typ_el, c_tensor, rho, 
     *    i_bnd_cdns, itermax, tol,
c    *    plot_modes, 
c    *    cmplx_max, real_max, int_max, 
     *    type_nod, 
     *    supplied_geo_flag, debug, show_mem_est, 
c     Inputs and outputs
     *    table_nod, type_el, x_arr,
c     Outputs
     *    v_nu_out, sol1, mode_pol, errco, emsg)

c         q_ac :   acoustic wave number (q_ac)
c         num_modes:  desired number of solved acoustic modes
c         n_msh_pts:  number of nodes in mesh
c         n_msh_el:   number of (triang) elements in mesh
c         n_type_el:  number of types of material 
c         type_nod:    ?? 
c         table_nod:    
c         type_el:
c         x_arr
c         v_nu_out:  eigen frequencies nu=omega/(2pi) for each mode
c         sol1:
c         mode_pol:  

C***********************************************************************
C
C  Program:
C     FEM solver of Acoustic waveguide problems.
C     This subroutine is compiled by f2py & called in mode_calcs.py
C
C  Authors:
C    Bjorn Sturmberg & Kokou B. Dossou
C
C***********************************************************************
C
      implicit none
C  Local parameters:
C       ! Propagation constant
      complex*16 q_ac 
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used
C     integer*8  plot_modes
      integer*8 errco
      character*2048 emsg

      integer :: stat=0
C       !  (int_max)
      integer*8, dimension(:), allocatable :: a   
C       !  (cmplx_max)
      complex*16, dimension(:), allocatable :: b   
C       !  (real_max)
      double precision, dimension(:), allocatable :: c   
      integer*8 supplied_geo_flag, symmetry_flag

C  Declare the pointers of the integer super-vector
      integer*8 ip_eq
      integer*8 ip_visite

C  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_mat2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_trav, jp_vp, jp_rhs
      integer*8 jp_eigenum_modes_tmp, jp_eigen_pol

c     Declare the pointers of the real super-vector
      integer*8 kp_mat1_re, kp_mat1_im
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im

c     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 ip_work, ip_work_sort, ip_work_sort2
      integer*8 nonz, nonz_max, max_row_len

      integer*8 n_typ_el
      complex*16 c_tensor(6,6,n_typ_el)
      complex*16 rho(n_typ_el)

      integer*8 i, j, ip
      integer*8 nodes_per_el, ui, debug, show_mem_est, namelength
      integer*8 n_msh_el, n_msh_pts, i_bnd_cdns, neq

C     ! Number of nodes per element
      parameter(nodes_per_el=6)
      integer*8 type_nod(n_msh_pts), type_el(n_msh_el) 
      integer*8 table_nod(nodes_per_el, n_msh_el)

      double precision pi
      double precision lat_vecs(2,2)
      double precision lx, ly, d_in_m

      complex*16 shift_nu, shift_omsq
      integer*8  i_base
      complex*16 ii

C  Variable used by valpr
      integer*8 ltrav, n_conv
      complex*16 z_beta, z_tmp, z_tmp0
      integer*8, dimension(:), allocatable :: iindex
c     variable used by UMFPACK

      double precision time1, time2
      double precision stime1, stime2
      character*(8) start_date, end_date
      character*(10) start_time, end_time

C  Variable used by valpr
      integer*8 num_modes, nvect, itermax
      double precision tol

C  Names and Controls
      character mesh_file*1000, gmsh_file*1000, log_file*1000
      character gmsh_file_pos*1000, dir_name*1000
      character*1000 tchar

c     new breed of variables to prise out of a, b and c
      double precision x_arr(2,n_msh_pts)
C       complex*16, target :: sol1(3,nodes_per_el+7,num_modes,n_msh_el)
      complex*16, target :: sol1(3,nodes_per_el,num_modes,n_msh_el)
      complex*16, target :: v_nu_out(num_modes)
      complex*16 mode_pol(4,num_modes)


Cf2py intent(in) q_ac, num_modes
Cf2py intent(in) debug, mesh_file, n_msh_pts, n_msh_el
Cf2py intent(in) d_in_m, shift
Cf2py intent(in) i_bnd_cdns, itermax, tol
Cf2py intent(in) plot_modes, c_tensor, rho
Cf2py intent(in) cmplx_max, real_max, int_max
Cf2py intent(in) n_typ_el, supplied_geo_flag
Cf2py intent(in) type_nod, table_nod, type_el, x_arr, symmetry_flag

C  Note: the dependent variables must be listed AFTER the
C  independent variables that they depend on in the function call!
Cf2py depend(c_tensor) n_typ_el
Cf2py depend(rho) n_typ_el
Cf2py depend(type_nod) n_msh_pts
Cf2py depend(table_nod) nodes_per_el, n_msh_el
Cf2py depend(type_el) n_msh_el
Cf2py depend(x_arr) n_msh_pts

Cf2py intent(out) v_nu_out
Cf2py intent(out) sol1, mode_pol, table_nod, type_el, x_arr

Cf2py intent(out) errco
Cf2py intent(out) emsg

C
CCCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCC
C
C     Set parameter for the super-vectors of integer and real numbers
C
C       !ui = Unite dImpression
      ui = 6     
C      nodes_per_el = 6 ! Number of nodes per element
      pi = 3.141592653589793d0
      ii = cmplx(0.0d0, 1.0d0, 8)

C       nvect = 2*num_modes + num_modes/2 +3
      nvect = 3*num_modes + 3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      errco= 0
      emsg = ""
 
c     Declare work space arrays

      call array_size(n_msh_pts, n_msh_el, num_modes, 
     *  int_max, cmplx_max, real_max, emsg, errco)
      RETONERROR(errco) 

      allocate(a(int_max), STAT=stat)
      call check_alloc(stat, int_max, "a", -1, errco, emsg)
      RETONERROR(errco) 

      allocate(b(cmplx_max), STAT=stat)
      call check_alloc(stat, cmplx_max, "b", -1, errco, emsg)
      RETONERROR(errco) 

      allocate(c(real_max), STAT=stat)
      call check_alloc(stat, real_max, "c", -1, errco, emsg)
      RETONERROR(errco) 

      allocate(iindex(num_modes), STAT=stat)
      call check_alloc(stat, num_modes, "iindex", -1, errco, emsg)
      RETONERROR(errco) 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     clean mesh_format
      namelength = len_trim(mesh_file)
      gmsh_file = mesh_file(1:namelength-5)//'.msh'
      gmsh_file_pos = mesh_file(1:namelength)
      log_file = mesh_file(1:namelength-5)//'-AC.log'
      if (debug .eq. 1) then
        write(*,*) "mesh_file = ", mesh_file
        write(*,*) "gmsh_file = ", gmsh_file
      endif

C       ! initial time  in unit = sec.
      call cpu_time(time1)  
      call date_and_time ( start_date, start_time )
C
C      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)
C      lx=1.0 ! Diameter of unit cell. Default, lx = 1.0.
C      ly=1.0 ! NOTE: currently requires ly=lx, ie rectangular unit cell
C     ToDo: sort out what's going on here - in EM lx=ly=1
      lx = d_in_m
      ly = d_in_m
      shift_omsq= (2*pi*shift_nu)**2

C
C####################  Start FEM PRE-PROCESSING  #######################
C

C       ! pointer to FEM connectivity table
      ip_visite = 1 
      ip_eq = ip_visite + n_msh_pts
      jp_x = 1
C
      if (supplied_geo_flag .eq. 0) then
        call geometry (n_msh_el, n_msh_pts, nodes_per_el, n_typ_el,
     *     lx, ly, type_nod, type_el, table_nod,
     *     x_arr, mesh_file, errco, emsg)
        if (errco .ne. 0) then
            return
        endif
      endif

      call lattice_vec (n_msh_pts, x_arr, lat_vecs, debug)

C       if (debug .eq. 1) then
C       open (unit=64, file="msh_check.txt",
C      *         status="unknown")
C         do i=1,n_msh_el
C           write(64,*) i, type_el(i)
C         enddo
C         write(64,*)
C         write(64,*)
C         write(64,*)
C         do i=1,n_msh_el
C           do j=1,nodes_per_el
C             write(64,*) i, j, table_nod(j,i)
C           enddo
C         enddo
C         write(64,*)
C         write(64,*)
C         write(64,*)
C         do j=1,nodes_per_el
C           write(64,*) j, type_nod(j)
C         enddo
C       close(63)
C       endif



      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: n_msh_pts, n_msh_el = ", 
     *    n_msh_pts, n_msh_el
      endif

c     Determine number of boundary conditions (neq) and 2D index array
c     a(ip_eq)
      call bound_cond_AC (i_bnd_cdns, n_msh_pts, neq, type_nod,
     *       a(ip_eq), debug)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Sparse matrix storage
      ip_col_ptr = ip_eq + 3*n_msh_pts

      call csr_max_length_AC (n_msh_el, n_msh_pts, neq, nodes_per_el,
     *  table_nod, a(ip_eq), a(ip_col_ptr), nonz_max)

      ip = ip_col_ptr + neq + 1
      if (ip .gt. int_max) then
         write(emsg,*) "py_calc_modes_AC: ip > int_max : ",
     *    ip, int_max,
     *    "py_calc_modes_AC: nonz_max = ", nonz_max,
     *    "py_calc_modes_AC: increase the size of int_max"
         errco = -3
         return 
      endif
c
      ip_row = ip_col_ptr + neq + 1

      call csr_length_AC (n_msh_el, n_msh_pts, neq, nodes_per_el, 
     *  table_nod, a(ip_eq), a(ip_row), a(ip_col_ptr), nonz_max,
     *  nonz, max_row_len, ip, int_max, debug)

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*n_msh_pts
      ip_work_sort2 = ip_work_sort + max_row_len

c     sorting csr ...
      call sort_csr (neq, nonz, max_row_len, a(ip_row),
     *  a(ip_col_ptr), a(ip_work_sort), a(ip_work),
     *  a(ip_work_sort2))

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: nonz_max = ", nonz_max
        write(ui,*) "py_calc_modes_AC: nonz = ", nonz
        write(ui,*) "py_calc_modes_AC: cmplx_max/nonz = ",
     *    dble(cmplx_max)/dble(nonz)
      endif

      int_used = ip_work_sort2 + max_row_len

      if (int_max .lt. int_used) then
        write(emsg,*)"The size of the integer supervector is too small",
     *   "integer super-vec: int_max  = ", int_max,
     *   "integer super-vec: int_used = ", int_used
        errco = -4
        return 
      endif

      jp_rhs = jp_x + 2*n_msh_pts
c     jp_rhs will also be used (in gmsh_post_process) to store a solution
      jp_mat2 = jp_rhs + max(neq, 3*n_msh_pts)
      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq
      jp_eigenum_modes_tmp = jp_resid+3*nodes_per_el*num_modes*n_msh_el
C       ! Eigenvectors
      jp_vschur = jp_eigenum_modes_tmp + num_modes + 1     
      jp_eigen_pol = jp_vschur + neq*nvect
      jp_trav = jp_eigen_pol + num_modes*4

      ltrav = 3*nvect*(nvect+2)
      jp_vp = jp_trav + ltrav
      cmplx_used = jp_vp + neq*num_modes

      if (cmplx_max .lt. cmplx_used)  then
        write(emsg,*)"The size of the complex supervector is too small",
     *   "complex super-vec: cmplx_max  = ", cmplx_max,
     *   "complex super-vec: cmplx_used = ", cmplx_used
         errco = -5
         return 
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
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

        write(emsg,*)"The size of the real supervector is too small",
     *   "2*nonz  = ", 2*nonz,
     *    "real super-vec: real_max  = ", real_max,
     *    "real super-vec: real_used = ", real_used

        errco = -6
        return 
      endif

c
c###############################################
c
c       ----------------------------------------------------------------
c       convert from 1-based to 0-based
c       ----------------------------------------------------------------
c
        do 60 j = 1, neq+1
            a(j+ip_col_ptr-1) = a(j+ip_col_ptr-1) - 1
60      continue
        do 70 j = 1, nonz
            a(j+ip_row-1) = a(j+ip_row-1) - 1
70      continue
c
c
c     The CSC iindexing, i.e., ip_col_ptr, is 1-based
c       (but valpr.f will change the CSC iindexing to 0-based iindexing)
      i_base = 0

C#####################  End FEM PRE-PROCESSING  #########################
C
      write(ui,*)
      write(ui,*) "-----------------------------------------------"
C       write(ui,*) " AC FEM, k_AC : ", real(q_ac), " 1/m"
C       write(ui,*) "-----------------------------------------------"
C       write(ui,*)

C       if (debug .eq. 1) then
C         write(ui,*) "py_calc_modes_AC: call to asmbly"
C       endif
      write(ui,*) "AC FEM: "
      write(ui,*) "      - assembling linear system"
C     Assemble the coefficient matrix K and M of the finite element equations

      call get_clocks(stime1, time1)

      call asmbly_AC (i_base, n_msh_el, n_msh_pts, neq, nodes_per_el,
     *  shift_omsq, q_ac, n_typ_el, rho, c_tensor,
     *  table_nod, type_el, a(ip_eq),
     *  x_arr, nonz, a(ip_row), a(ip_col_ptr),
     *  c(kp_mat1_re), c(kp_mat1_im), b(jp_mat2), a(ip_work), 
     *  symmetry_flag, debug)

      call get_clocks(stime2, time2)
      write(ui,'(A,F6.2,A)') '           cpu time  = ', (time2-time1), 
     *   ' secs.'
      write(ui,'(A,F6.2,A)') '           wall time = ', (stime2-stime1), 
     *   ' secs.'

C       if (debug .eq. 1) then
C         write(ui,*) "py_calc_modes_AC: call to valpr"
C       endif
      write(ui,*) "      - solving linear system"
      call get_clocks(stime1, time1)


      call valpr_64_AC (i_base, nvect, num_modes, neq, itermax, ltrav,
     *  tol, nonz, a(ip_row), a(ip_col_ptr), c(kp_mat1_re),
     *  c(kp_mat1_im), b(jp_mat2), 
     *  b(jp_vect1), b(jp_vect2), b(jp_workd), b(jp_resid), 
     *     b(jp_vschur), v_nu_out, b(jp_trav), b(jp_vp), 
     *  c(kp_rhs_re), c(kp_rhs_im), c(kp_lhs_re), c(kp_lhs_im), n_conv, 
     *  debug, show_mem_est, errco, emsg)


      if (errco .ne. 0) then
          return
      endif

      call get_clocks(stime2, time2)
      write(ui,'(A,F6.2,A)') '           cpu time  = ', (time2-time1), 
     *   ' secs.'
      write(ui,'(A,F6.2,A)') '           wall time = ', (stime2-stime1),
     *   ' secs.'


      if (n_conv .ne. num_modes) then
         write(emsg, '(A,I5,I5)') 
     *    "py_calc_modes_AC: convergence problem " //
     *    "in valpr_64: n_conv != num_modes  ", n_conv, num_modes
         errco = -7
         return 
      endif
C
      do i=1,num_modes
        z_tmp0 = v_nu_out(i)
        z_tmp = 1.0d0/z_tmp0+shift_omsq
        z_beta = sqrt(z_tmp) / (2.0d0 * pi)
C       Frequency (z_beta) should always be positive.
        if (dble(z_beta) .lt. 0) z_beta = -z_beta
        v_nu_out(i) = z_beta
      enddo
c
      call z_indexx_AC (num_modes, v_nu_out, iindex)
C
C       The eigenvectors will be stored in the array sol1
C       The eigenum_modesues and eigenvectors will be renumbered
C                 using the permutation vector iindex
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to array_sol"
      endif
        call array_sol_AC (num_modes, n_msh_el, n_msh_pts, neq, 
     *   nodes_per_el, iindex, table_nod, type_el, a(ip_eq), x_arr, 
     *   v_nu_out,  b(jp_eigenum_modes_tmp), mode_pol, b(jp_vp), sol1)

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: array_sol returns call"
      endif
C
      if(debug .eq. 1) then
        write(ui,*) 'iindex = ', (iindex(i), i=1,num_modes)
      endif
      if(debug .eq. 1) then
        write(ui,*)
C         write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
C         write(ui,*) "sqrt(shift_omsq)/(2*pi) = ", sqrt(omsq) / (2.0d0 * pi)
        do i=1,num_modes
          write(ui,"(i4,2(g22.14),2(g18.10))") i, v_nu_out(i)
        enddo
      endif
      
CC    Save Original solution
C      if (plot_modes .eq. 1) then
C        dir_name = "AC_fields"
CC        call write_sol_AC (num_modes, n_msh_el, nodes_per_el, lambda,
CC      *       v_nu_out, sol1, mesh_file, dir_name)
CC        call write_param (lambda, n_msh_pts, n_msh_el, i_bnd_cdns,
CC    *       num_modes, nvect, itermax, tol, shift_omsq, lx, ly,
CC    *       mesh_file, n_conv, dir_name)
C        tchar = "AC_fields/All_plots_png_abs2_eE.geo"
C        open (unit=34,file=tchar)
C          do i=1,num_modes
C            call gmsh_post_process_AC (i, num_modes, n_msh_el, 
C     *         n_msh_pts, nodes_per_el, table_nod, type_el,
C     *         x_arr, v_nu_out, sol1, b(jp_rhs), a(ip_visite),
C     *         gmsh_file_pos, dir_name, d_in_m, debug)
C          enddo
C        close (unit=34)
C      endif
CC
C
C#########################  End Calculations  ###########################
C
      call date_and_time ( end_date, end_time )
      call cpu_time(time2)
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) 'Total CPU time (sec.)  = ', (time2-time1)
C
        open (unit=26,file=log_file)
        write(26,*)
        write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
        write(26,*) "Start date and time   = ", start_date,
     *    " ; ", start_time
        write(26,*) "End date and time     = ", end_date,
     *    " ; ", end_time
        write(26,*) "Total CPU time (sec.) = ",  (time2-time1)
        write(26,*)
        write(26,*) "q_ac = ", q_ac
        write(26,*) "shift_omsq= ", shift_omsq
        write(26,*)
        write(26,*) "n_msh_pts, n_msh_el, nodes_per_el  = ", n_msh_pts, 
     *               n_msh_el, nodes_per_el
        write(26,*) "neq, i_bnd_cdns = ", neq, i_bnd_cdns
        write(26,*) " lat_vecs:  = "
        write(26,"(2(f18.10))") lat_vecs
        write(26,*) "mesh_file = ", mesh_file
        write(26,*) "gmsh_file = ", gmsh_file
        write(26,*) "log_file  = ", log_file
        close(26)
C
        write(ui,*) "   .      .      ."
        write(ui,*) "   .      .      ."
        write(ui,*) "   .      . (d=",d_in_m,")"
        write(ui,*) "  and   we're   done!"
      endif

      write(ui,*) "-----------------------------------------------"
      write(ui,*)
C
      deallocate(a,b,c,iindex)

      end subroutine calc_AC_modes

