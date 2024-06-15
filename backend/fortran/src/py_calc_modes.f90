#include "numbat_decl.h"


!  lambda - free space wavelength in m
!  n_modes - desired number of eigenvectors
!  n_msh_pts - number of FEM mesh points
!  n_msh_el  - number of FEM (triang) elements
!  n_typ_el  - number of types of elements (and therefore elements)
!  v_refindex_n - array of effective index of materials
!  bloch_vec - in-plane k-vector (normally tiny just to avoid degeneracies)
!  shift_ksqr   - k_est^2 = n^2 k_0^2  : estimate of eigenvalue k^2
!  bnd_cnd_i - bnd conditions (Dirichlet = 0, Neumann = 1, Periodic = 2)
!  v_eigs_beta  - array of eigenvalues kz
!  sol1   - 4-dim array of solutions [field comp, node of element (1..13)?!, eigvalue, element number] (strange ordering)
!  mode_pol  - unknown - never used in python
!  table_nod - 2D array [node_on_elt-1..6][n_msh_el] giving the mesh point mp of each node
!  Points where type_el[mp] is not the same for all 6 nodes must be interface points
!  type_el   - n_msh_el array: material index for each element
!  type_nod  - is boundary node?
!  mesh_xy  - (2 , n_msh_pts)  x,y coords?
!  ls_material  - (1, nodes_per_el+7, n_msh_el)

subroutine calc_EM_modes( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &    ! inputs
   E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_typ_el, v_refindex_n, & ! inputs
   v_eigs_beta, sol1, mode_pol, table_nod, type_el, type_nod, mesh_xy, ls_material, errco, emsg)  ! output params

   use numbatmod

   use nbinterfaces


   integer*8, parameter :: nodes_per_el = 6

   double precision lambda, dimscale_in_m, bloch_vec(2)

   integer*8 n_modes, n_typ_el
   integer*8 n_msh_el, n_msh_pts,  bdy_cdn, itermax
   integer*8 neq, debug
   complex*16 shift_ksqr

   integer*8 type_el(n_msh_el), type_nod(n_msh_pts)
   integer*8 table_nod(nodes_per_el, n_msh_el)
   integer*8 E_H_field

   complex*16 ls_material(1,nodes_per_el+7,n_msh_el)
   double precision mesh_xy(2,n_msh_pts)
   complex*16 mode_pol(4,n_modes)


   complex*16  v_refindex_n(n_typ_el)
   integer errco
   character(len=EMSG_LENGTH) emsg




   integer*8 int_max, cmplx_max, int_used, cmplx_used
   integer*8 real_max, real_used

   integer*8, dimension(:), allocatable :: a_iwork
   complex*16, dimension(:), allocatable :: b_zwork
   double precision, dimension(:), allocatable :: c_dwork
   double precision, dimension(:,:), allocatable :: d_dwork
   integer*8, dimension(:), allocatable :: iindex
   complex*16, dimension(:,:), allocatable :: overlap_L
!
!  Declare the pointers of the integer super-vector
   integer*8 ip_table_E, ip_table_N_E_F, ip_visited
   integer*8 ip_type_N_E_F, ip_eq
   integer*8 ip_period_N, ip_nperiod_N
   integer*8 ip_period_N_E_F, ip_nperiod_N_E_F

!  Declare the pointers of the real super-vector
   integer*8 jp_x_N_E_F

   integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
   integer*8 jp_trav, jp_vp
   complex*16 pp(n_typ_el), qq(n_typ_el)
   complex*16 eps_eff(n_typ_el)
!





   integer*8 n_msh_pts_p3, ui

!  Variable used by valpr
   integer*8 nvect, ltrav
   integer*8 n_conv, i_base
   double precision ls_data(10)

   integer*8 n_core(2)  ! index of highest epsilon material, seems funky
   complex*16 z_beta, z_tmp, z_tmp0
   integer*8 n_edge, n_face, n_ddl, n_ddl_max, n_k

!  variable used by UMFPACK
   double precision control (20), info_umf (90)
   integer*8 numeric



   integer*8 i
   integer*8 ival, iel, inod


   double precision freq, tol
   double precision k_0, dim_x, dim_y
   double precision  bloch_vec_k(2)


!  Timing variables
   double precision time1, time2
   double precision stime1, stime2
   double precision time1_fact, time2_fact
   double precision time1_postp
   double precision stime1_postp
   double precision time1_arpack, time2_arpack
   double precision time1_J, time2_J
   double precision stime1_J, stime2_J

   character*(10) start_time, end_time

   !Names and Controls
   character(len=FNAME_LENGTH)  mesh_file, gmsh_file, log_file, gmsh_file_pos, overlap_file

   character msg*20
   integer*8 namelength
   integer*8 pair_warning, homogeneous_check

!  Declare the pointers of the real super-vector
   integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im
   integer*8 kp_mat1_re, kp_mat1_im

!  Declare the pointers of for sparse matrix storage
   integer*8 ip_col_ptr, ip_row
   integer*8 jp_mat2
   integer*8 ip_work, ip_work_sort, ip_work_sort2
   integer*8 nonz, nonz_max, max_row_len

   integer*8 ip

!  new breed of variables to prise out of a_iwork, b_zwork and c_dwork

   complex*16, target :: sol1(3,nodes_per_el+7,n_modes,n_msh_el)
   complex*16, target :: sol2(3,nodes_per_el+7,n_modes,n_msh_el)
   complex*16, pointer :: sol(:,:,:,:)


   complex*16, target :: v_eigs_beta(n_modes), beta2(n_modes)
   complex*16, pointer :: beta(:)



   integer *8:: ilo, ihi

   integer :: is_em



!f2py intent(in) lambda, n_modes
!f2py intent(in) debug, mesh_file, n_msh_pts, n_msh_el
!f2py intent(in) v_refindex_n, bloch_vec, dimscale_in_m, shift_ksqr
!f2py intent(in) E_H_field, bdy_cdn, itermax
!f2py intent(in) plot_modes, plot_real, plot_imag, plot_abs
!f2py intent(in) cmplx_max, real_max, int_max, n_typ_el

!f2py depend(v_refindex_n) n_typ_el

!f2py intent(out) v_eigs_beta, type_nod, ls_material
!f2py intent(out) sol1, mode_pol, table_nod, type_el, mesh_xy

!f2py intent(out) errco
!f2py intent(out) emsg



   tol = 0.d0

   errco = 0
   emsg = ""

!  Declare work space arrays

   is_em = 1

   !TODO: make a new dwork for x_N_E_F
   call prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, int_max, cmplx_max, real_max, &
      a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L,  errco, emsg)
   RETONERROR(errco)



!CCCCCCCCCCCCCCCC POST F2PY CCCCCCCCCCCCCCCCCCCCCCCCC

!  clean mesh_format
!  TODO: Drop these debug files?
   namelength = len_trim(mesh_file)
   gmsh_file = mesh_file(1:namelength-5)//'.msh'
   gmsh_file_pos = mesh_file(1:namelength)
   log_file = mesh_file(1:namelength-5)//'.log'
   if (debug .eq. 1) then
      write(*,*) "mesh_file = ", mesh_file
      write(*,*) "gmsh_file = ", gmsh_file
   endif

!  Calculate effective permittivity for each type of material
   !do i_32 = 1, int(n_typ_el)
   !   eps_eff(i_32) = v_refindex_n(i_32)**2
   !end do
   eps_eff = v_refindex_n**2

!  !ui = Unite dImpression
   ui = 6

!  nsym = 1 ! nsym = 0 => symmetric or hermitian matrices
!

   nvect = 2*n_modes + n_modes/2 +3

!CCCCCCCCCCCCCCCC END POST F2PY CCCCCCCCCCCCCCCCCCCCC

!  ! initial time  in unit = sec.
   call get_clocks(stime1, time1)
!
!####################  Start FEM PRE-PROCESSING  ########################
!

   dim_x = dimscale_in_m
   dim_y = dimscale_in_m

   ! Fill:  mesh_xy, type_nod, type_el, table_nod
   call geometry (n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, dim_x, dim_y, mesh_file, &
      mesh_xy, type_nod, type_el, table_nod,  errco, emsg)
   RETONERROR(errco)


   ! Storage locations in sequence
   !  - table_edge_face = a_iwork(ip_table_N_E_F),   shape: 14 x n_msh_el
   !  - visited         = a_iwork(ip_visited),       shape: n_ddl_max ? = n_msh_pts + n_msh_el
   !  - table_edges     = a_iwork(ip_table_E)        shape: 4 x n_msh_pts
   !
   !   visited is used as workspace. has no meaning between functions
   !
   !   V = number of vertices
   !   E = number of edges
   !   F = number of faces
   !   C = number of cells (3D, tetrahedron)
   !
   !  From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
   !  n_msh_pts = (number of vertices) + (number of mid-edge point) = V + E;
   !
   !

   ! Fills:  table_edge_face[1,:]
   ip_table_N_E_F = 1
   call list_face (n_msh_el, a_iwork(ip_table_N_E_F))

   !  n_ddl_max = max(N_Vertices) + max(N_Edge) + max(N_Face)
   !  For P2 FEM n_msh_pts=N_Vertices+N_Edge
   !  note: each element has 1 face, 3 edges and 10 P3 nodes
   !        so table_N_E_F = table_edge_face has dimensions 14 x n_msh_el

   ! each element is a face
   n_face = n_msh_el

   n_ddl_max = n_msh_pts + n_face

   ip_visited =  ip_table_N_E_F  + 14*n_msh_el
   ip_table_E = ip_visited + n_ddl_max

   ! Fills: n_edge, table_edge[1..4,:], table_edge_face[2:4,:], visited[1:n_msh_pts]
   call list_edge (n_msh_el, n_msh_pts, nodes_per_el, n_edge, type_nod, table_nod, &
      a_iwork(ip_table_E), a_iwork(ip_table_N_E_F), a_iwork(ip_visited))

   ! Fills: remainder of table_edge_face[5:,:], visited[1:n_msh_pts]
   call list_node_P3 (n_msh_el, n_msh_pts, nodes_per_el, n_edge, n_msh_pts_p3, table_nod, &
      a_iwork(ip_table_N_E_F), a_iwork(ip_visited))

   ! TODO: what is signif of this quanitty?
   n_ddl = n_edge + n_face + n_msh_pts_p3


   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
      write(ui,*) "py_calc_modes.f: n_msh_pts_p3 = ", n_msh_pts_p3
      write(ui,*) "py_calc_modes.f: n_vertex, n_edge, n_face,", " n_msh_el = ", &
         (n_msh_pts - n_edge), n_edge, n_face, n_msh_el
      write(ui,*) "py_calc_modes.f: 2D case of the Euler &
      & characteristic: V-E+F=1-(number of holes)"
      write(ui,*) "py_calc_modes.f: Euler characteristic: V - E + F &
      &= ", (n_msh_pts - n_edge) - n_edge + n_face
   endif
!C
!C-----------
!C
!C  overwriting pointers ip_row_ptr, ..., ip_adjncy
!
   ip_type_N_E_F = ip_table_E + 4*n_edge   ! not sure why 4* n_edge, not 4*n_msh_pts?
!

   ! TODO:
   !  ip is an index into an a_iwork, make this clearer!
   !  jp is an index into an b_zwork
   !  kp is an index into an c_dwork

   ! Offsets into the b_zwork workspace
   jp_x_N_E_F = 1


   ! Offsets into the a_iwork workspace
   ip_period_N = ip_type_N_E_F + 2*n_ddl
   ip_nperiod_N = ip_period_N + n_msh_pts
   ip_period_N_E_F = ip_nperiod_N + n_msh_pts
   ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
   ip_eq = ip_nperiod_N_E_F + n_ddl


   ! Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   ! Should be using c_dwork for x_E_F ?
   call type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, type_nod, table_nod, &
      a_iwork(ip_table_N_E_F), a_iwork(ip_visited), a_iwork(ip_type_N_E_F), mesh_xy, &
      !b_zwork(jp_x_N_E_F) &
   d_dwork &
      )


   ! Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   call get_coord_p3 (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, table_nod, type_nod, &
      a_iwork(ip_table_N_E_F), a_iwork(ip_type_N_E_F), mesh_xy, &
      !b_zwork(jp_x_N_E_F), &
   d_dwork, &
      a_iwork(ip_visited))




   ! TODO: the b_zwork should actually be the d_dwork containing x_N_E_F, but only matters for periodic
   call set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el, mesh_xy, nodes_per_el, &
      type_nod, table_nod, n_ddl, neq, ip_type_N_E_F, ip_eq, a_iwork, &
      !b_zwork, &
   d_dwork,  &
      int_max, debug)




!
!  Sparse matrix storage
!
   ip_col_ptr = ip_eq + 3*n_ddl

   call csr_max_length (n_msh_el, n_ddl, neq, nodes_per_el, a_iwork(ip_table_N_E_F), &
      a_iwork(ip_eq), a_iwork(ip_col_ptr), nonz_max)
!
!  ip = ip_col_ptr + neq + 1 + nonz_max
   ip = ip_col_ptr + neq + 1
   if (ip .gt. int_max) then
      write(emsg,*) "py_calc_modes.f: ip > int_max : ", ip, int_max, "py_calc_modes.f: nonz_max = ", &
         nonz_max, "py_calc_modes.f: increase the size of int_max"
      errco = -11
      return
   endif
!
   ip_row = ip_col_ptr + neq + 1

   call csr_length (n_msh_el, n_ddl, neq, nodes_per_el, a_iwork(ip_table_N_E_F), a_iwork(ip_eq), a_iwork(ip_row), &
      a_iwork(ip_col_ptr), nonz_max, nonz, max_row_len, ip, int_max, debug)

   ip_work = ip_row + nonz
   ip_work_sort = ip_work + 3*n_ddl
   ip_work_sort2 = ip_work_sort + max_row_len

!  sorting csr ...
   call sort_csr (neq, nonz, max_row_len, a_iwork(ip_row), a_iwork(ip_col_ptr), a_iwork(ip_work_sort), a_iwork(ip_work), &
      a_iwork(ip_work_sort2))

   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: nonz_max = ", nonz_max
      write(ui,*) "py_calc_modes.f: nonz = ", nonz
      write(ui,*) "py_calc_modes.f: cmplx_max/nonz = ", dble(cmplx_max)/dble(nonz)
   endif

   int_used = ip_work_sort2 + max_row_len

   if (int_max .lt. int_used) then
      write(emsg,*)'The size of the integer supervector is too small', 'integer super-vec: int_max  = ', &
         int_max, 'integer super-vec: int_used = ', int_used
      errco = -12
      return
   endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!
!C
   jp_mat2 = jp_x_N_E_F + 3*n_ddl

   jp_vect1 = jp_mat2 + nonz
   jp_vect2 = jp_vect1 + neq
   jp_workd = jp_vect2 + neq
   jp_resid = jp_workd + 3*neq

!  ! Eigenvectors
   jp_vschur = jp_resid + neq
   jp_trav = jp_vschur + neq*nvect

   ltrav = 3*nvect*(nvect+2)
   jp_vp = jp_trav + ltrav

   cmplx_used = jp_vp + neq*n_modes
!
   if (cmplx_max .lt. cmplx_used)  then
      write(emsg,*)'The size of the complex supervector is too small', 'complex super-vec: int_max  = ', &
         cmplx_max, 'complex super-vec: int_used = ', cmplx_used
      errco = -13
      return
   endif
!
   kp_rhs_re = 1
   kp_rhs_im = kp_rhs_re + neq
   kp_lhs_re = kp_rhs_im + neq
   kp_lhs_im = kp_lhs_re + neq
   kp_mat1_re = kp_lhs_im + neq
   kp_mat1_im = kp_mat1_re + nonz
   real_used = kp_mat1_im + nonz

   if (real_max .lt. real_used) then
      write(emsg,*) 'The size of the real supervector is too small', '2*nonz  = ', 2*nonz, &
         'real super-vec: real_max  = ', real_max, 'real super-vec: real_used = ', real_used
      errco = -14
      return
   endif
!
!
!###############################################
!
!  ----------------------------------------------------------------
!  convert from 1-based to 0-based
!  ----------------------------------------------------------------
!
! do j = 1, neq+1
! a_iwork(j+ip_col_ptr-1) = a_iwork(j+ip_col_ptr-1) - 1
! end do
! do  j = 1, nonz
! a_iwork(j+ip_row-1) = a_iwork(j+ip_row-1) - 1
! end do
! !

   ilo = ip_col_ptr-1 + 1
   ihi = ip_col_ptr-1 + neq + 1
   a_iwork(ilo:ihi) = a_iwork(ilo:ihi) - 1

   ilo = ip_row-1 + 1
   ihi = ip_row-1 + nonz
   a_iwork(ilo:ihi) = a_iwork(ilo:ihi) - 1




!  The CSC indexing, i.e., ip_col_ptr, is 1-based
!  (but valpr.f will change the CSC indexing to 0-based indexing)
   i_base = 0
!
!
!#####################  End FEM PRE-PROCESSING  #########################
!
!
!
   write(ui,*)
   write(ui,*) "-----------------------------------------------"
!  write(ui,*) " EM FEM, lambda : ", lambda*1.0d9, "nm"
!  write(ui,*) "-----------------------------------------------"
!  write(ui,*)
!  C
   freq = 1.0d0/lambda
   k_0 = 2.0d0*D_PI*freq
!
!  Index number of the core materials (material with highest Re(eps_eff))
   if(dble(eps_eff(1)) .gt. dble(eps_eff(2))) then
      n_core(1) = 1
   else
      n_core(1) = 2
   endif
   n_core(2) = n_core(1)

!  Check that the layer is not in fact homogeneous
   homogeneous_check = 0
   do i=1,n_typ_el-1
      if(dble(eps_eff(i)) .ne. dble(eps_eff(i+1))) then
         homogeneous_check = 1
      elseif(dimag(eps_eff(i)) .ne. dimag(eps_eff(i+1))) then
         homogeneous_check = 1
      endif
   enddo

   if(homogeneous_check .eq. 0) then
      write(emsg,*) "py_calc_modes_1d.f: ", "FEM routine cannot handle homogeneous layers.", &
         "Define layer as object.ThinFilm"
      errco = -17
      return
   endif
!  Parameter for shift-and-invert method - now given as input from python
!  shift_ksqr = 1.01d0*Dble(v_refindex_n(n_core(1)))**2 * k_0**2
!  * - bloch_vec(1)**2 - bloch_vec(2)**2


   if(debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: n_core = ", n_core
      if(E_H_field .eq. FEM_FORMULATION_E) then
         write(ui,*) "py_calc_modes.f: E-Field formulation"
      else
         write(ui,*) "py_calc_modes.f: H-Field formulation"
      endif
   endif
!
!
   if(E_H_field .eq. FEM_FORMULATION_E) then
      do i=1,n_typ_el
         qq(i) = eps_eff(i)*k_0**2
         pp(i) = 1.0d0
      enddo
   elseif(E_H_field .eq. FEM_FORMULATION_H) then
      do i=1,n_typ_el
         qq(i) = k_0**2
         pp(i) = 1.0d0/eps_eff(i)
      enddo
   else
!  Can't happen
      write(ui,*) "py_calc_modes.f: action indef. avec E_H_field = ", E_H_field
      write(ui,*) "Aborting..."
      errco = -18
      return
   endif
!
!CCCCCCCCCCCCCCCCCCC  Loop over Adjoint and Prime  CCCCCCCCCCCCCCCCCCCCCC
!
   do n_k = 1,2
!
      if (n_k .eq. 1) then
         sol => sol1
         beta => v_eigs_beta
         bloch_vec_k = bloch_vec
         msg = "adjoint solution"
      else
         sol => sol2
         beta => beta2
         bloch_vec_k = -bloch_vec
         msg = "prime solution"
      endif
!
!  Assemble the coefficient matrix A and the right-hand side F of the
!  finite element equations
!  if (debug .eq. 1) then
!  write(ui,*) "py_calc_modes.f: Asmbly: call to asmbly"
!  endif
      write(ui,*) "EM FEM: "
      write(ui,'(A,A)') "   - assembling linear system for ", msg

      call get_clocks(stime1, time1)

      call asmbly ( bdy_cdn, i_base, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, shift_ksqr, &
         bloch_vec_k, n_typ_el, pp, qq, table_nod, a_iwork(ip_table_N_E_F), type_el, a_iwork(ip_eq), &
         a_iwork(ip_period_N), a_iwork(ip_period_N_E_F), mesh_xy, &
         !b_zwork(jp_x_N_E_F), &
      d_dwork, &
         nonz, a_iwork(ip_row), &
         a_iwork(ip_col_ptr), c_dwork(kp_mat1_re), c_dwork(kp_mat1_im), b_zwork(jp_mat2), a_iwork(ip_work))

      call get_clocks(stime2, time2)
      write(ui,'(A,F6.2,A)') '  cpu time  = ', (time2-time1), ' secs.'
      write(ui,'(A,F6.2,A)') '  wall time = ', (stime2-stime1), ' secs.'
!
!  factorization of the globale matrice
!  -----------------------------------
!
      if (debug .eq. 1) then
         write(ui,*) "py_calc_modes.f: Adjoint(1) / Prime(2)", n_k
!  write(ui,*) "py_calc_modes.f: factorisation: call to znsy"
      endif
!
!  if (debug .eq. 1) then
!  write(ui,*) "py_calc_modes.f: call to valpr"
!  endif
      write(ui,*) "  - solving linear system"

      call get_clocks(stime1, time1)

      call valpr_64 (i_base, nvect, n_modes, neq, itermax, ltrav, tol, nonz, a_iwork(ip_row), a_iwork(ip_col_ptr), &
         c_dwork(kp_mat1_re), c_dwork(kp_mat1_im), b_zwork(jp_mat2), b_zwork(jp_vect1), b_zwork(jp_vect2), b_zwork(jp_workd), &
         b_zwork(jp_resid), b_zwork(jp_vschur), beta, b_zwork(jp_trav), b_zwork(jp_vp), c_dwork(kp_rhs_re), c_dwork(kp_rhs_im), &
         c_dwork(kp_lhs_re), c_dwork(kp_lhs_im), n_conv, ls_data, numeric, control, info_umf, debug, errco, emsg)
      call get_clocks(stime2, time2)
      write(ui,'(A,F6.2,A)') '  cpu time  = ', (time2-time1), ' secs.'
      write(ui,'(A,F6.2,A)') '  wall time = ', (stime2-stime1), ' secs.'
!
      if (errco .ne. 0) then
         return
      endif

      if (n_conv .ne. n_modes) then
         write(ui,*) "py_calc_modes.f: convergence problem in valpr_64"
!  write(ui,*) "You should probably increase resolution of mesh!"
         write(ui,*) "py_calc_modes.f: n_conv != n_modes : ", n_conv, n_modes
!  write(ui,*) "n_core(1), v_refindex_n(n_core(1)) = ",
!  * n_core(1), v_refindex_n(n_core(1))
         write(ui,*) "py_calc_modes.f: Aborting..."
         errco = -19
         return
      endif
!
      time1_fact = ls_data(1)
      time2_fact = ls_data(2)
!
      time1_arpack = ls_data(3)
      time2_arpack = ls_data(4)
!
      do i=1,n_modes
         z_tmp0 = beta(i)
         z_tmp = 1.0d0/z_tmp0+shift_ksqr
         z_beta = sqrt(z_tmp)
!  Mode classification - we want the forward propagating mode
         if (abs(imag(z_beta)/z_beta) .lt. 1.0d-8) then
!  re(z_beta) > 0 for forward propagating mode
            if (dble(z_beta) .lt. 0) z_beta = -z_beta
         else
!  im(z_beta) > 0 for forward decaying evanescent mode
            if (imag(z_beta) .lt. 0) z_beta = -z_beta
         endif
!  !  Effective iindex
!  z_beta = sqrt(z_tmp)/k_0
         beta(i) = z_beta
      enddo
!
      call get_clocks(stime1_postp, time1_postp)
!

!  order beta by magnitudes and store in iindex
      call z_indexx (n_modes, beta, iindex)
!
!  The eigenvectors will be stored in the array sol
!  The eigenvalues and eigenvectors are renumbered
!  using the permutation vector iindex
      call array_sol ( bdy_cdn, n_modes, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, n_core, &
         bloch_vec_k, iindex, table_nod, a_iwork(ip_table_N_E_F), type_el, a_iwork(ip_eq), a_iwork(ip_period_N), &
         a_iwork(ip_period_N_E_F), mesh_xy, &
      !b_zwork(jp_x_N_E_F),
         d_dwork, &
         beta, mode_pol, b_zwork(jp_vp), sol)
!
      if(debug .eq. 1) then
         write(ui,*) 'iindex = ', (iindex(i), i=1,n_modes)
      endif
      if(debug .eq. 1) then
         write(ui,*)
         write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
         write(ui,*) (bloch_vec_k(i)/(2.0d0*D_PI),i=1,2)
         write(ui,*) "sqrt(shift_ksqr) = ", sqrt(shift_ksqr)
         write(ui,*) "n_modess = "
         do i=1,n_modes
            write(ui,"(i4,2(g22.14),2(g18.10))") i, beta(i)
         enddo
      endif
!
!  Calculate energy in each medium (typ_el)
      if (n_k .eq. 2) then
         call mode_energy (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_core, table_nod, type_el, n_typ_el, eps_eff,&
            mesh_xy, sol, beta, mode_pol)
      endif
!
   enddo
!
!CCCCCCCCCCCCCCCCCCCCCCC  End Prime, Adjoint Loop  CCCCCCCCCCCCCCCCCCCCCC
!
!  Orthogonal integral
   pair_warning = 0
   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: Field product"
   endif
   overlap_file = "Orthogonal.txt"
   call get_clocks(stime1_J, time1_J)
   call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, pp, table_nod, &
      type_el, mesh_xy, v_eigs_beta, beta2, &
      sol1, sol2, overlap_L, overlap_file, debug, pair_warning, k_0)

   if (pair_warning .ne. 0 .and. n_modes .le. 20) then
      write(ui,*) "py_calc_modes.f: Warning found 1 BM of cmplx conj"
      write(ui,*) "pair, increase num_BMs to include the other."
   endif

   call get_clocks(stime2_J, time2_J)
   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: CPU time for orthogonal :", (time2_J-time1_J)
   endif

!  The z-component must be multiplied by -ii*beta in order to
!  get the physical, un-normalised z-component
!  (see Eq. (25) of the JOSAA 2012 paper)
   do ival=1,n_modes
      do iel=1,n_msh_el
         do inod=1,nodes_per_el+7
            sol1(3,inod,ival,iel) = C_IM_ONE * beta(ival) * sol1(3,inod,ival,iel)
         enddo
      enddo
   enddo

   call array_material_EM (n_msh_el, n_typ_el, v_refindex_n, type_el, ls_material)

!
!  Save Original solution
!  if (plot_modes .eq. 1) then
!  dir_name = "Bloch_fields"
!  q_average = 0
!C  call write_sol (n_modes, n_msh_el, nodes_per_el, E_H_field, lambda,
!C  * v_eigs_beta, sol1, mesh_file, dir_name)
!C  call write_param (E_H_field, lambda, n_msh_pts, n_msh_el, bdy_cdn,
!C  * n_modes, nvect, itermax, tol, shift_ksqr, dim_x, dim_y,
!C  * mesh_file, mesh_format, n_conv, n_typ_el, eps_eff,
!C  * bloch_vec, dir_name)
!  tchar = "Bloch_fields/PDF/All_plots_pdf.geo"
!  open (unit=34,file=tchar)
!  do i=1,n_modes
!  call gmsh_post_process (i, E_H_field, n_modes, n_msh_el,
!  *  n_msh_pts, nodes_per_el, table_nod, type_el, n_typ_el,
!  *  v_refindex_n, mesh_xy, v_eigs_beta, sol1,
!  *  a_iwork(ip_visited), gmsh_file_pos, dir_name,
!  *  q_average, plot_real, plot_imag, plot_abs)
!  enddo
!  close (unit=34)
!  endif
!C

!  Normalisation
   if(debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: Field  Normalisation"
   endif
   call get_clocks(stime1_J, time1_J)
   call normalisation (n_modes, n_msh_el, nodes_per_el, sol1, sol2, overlap_L)
   call get_clocks(stime2_J, time2_J)
   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: CPU time for normalisation :", (time2_J-time1_J)
   endif
!
!  Orthonormal integral
   if (debug .eq. 1) then
      write(ui,*) "py_calc_modes.f: Product of normalised field"
      overlap_file = "Orthogonal_n.txt"
      call get_clocks(stime1_J, time1_J)
      call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, pp, table_nod, &
         type_el, mesh_xy, v_eigs_beta, beta2, sol1, sol2, overlap_L, overlap_file, debug, &
         pair_warning, k_0)
      call get_clocks(stime2_J, time2_J)
      write(ui,*) "py_calc_modes.f: CPU time for orthogonal :", (time2_J-time1_J)
   endif
!
!#########################  End Calculations  ###########################
!
   call get_clocks(stime2, time2)
!
   if (debug .eq. 1) then
      write(ui,*)
      write(ui,*) 'Total CPU time (sec.)  = ', (time2-time1)

      open (unit=26,file=log_file)
      write(26,*)
      write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
      write(26,*) "Start time   = ", start_time
      write(26,*) "End time  = ", end_time
      write(26,*) "Total CPU time (sec.) = ", (time2-time1)
      write(26,*) "LU factorisation : CPU time and % Total time = ", (time2_fact-time1_fact), &
         100*(time2_fact-time1_fact)/(time2-time1),"%"
      write(26,*) "ARPACK : CPU time and % Total time = ", (time2_arpack-time1_arpack), &
         100*(time2_arpack-time1_arpack)/(time2-time1),"%"
!  write(26,*) "Assembly : CPU time and % Total time = ",
!  *   (time2_asmbl-time1_asmbl),
!  *   100*(time2_asmbl-time1_asmbl)/(time2-time1),"%"
      write(26,*) "Post-processsing : CPU time and % Total time = ", (time2-time1_postp), &
         100*(time2-time1_postp)/(time2-time1),"%"
!  write(26,*) "Pre-Assembly : CPU time and % Total time = ",
!  *   (time1_asmbl-time1),
!  *   100*(time1_asmbl-time1)/(time2-time1),"%"
      write(26,*)
      write(26,*) "lambda  = ", lambda
      write(26,*) "n_msh_pts, n_msh_el, nodes_per_el  = ", n_msh_pts, n_msh_el, nodes_per_el
      write(26,*) "neq, bdy_cdn = ", neq, bdy_cdn
      if ( E_H_field .eq. FEM_FORMULATION_E) then
         write(26,*) "E_H_field   = ", E_H_field, " (E-Field formulation)"
      elseif ( E_H_field .eq. FEM_FORMULATION_H) then
         write(26,*) "E_H_field   = ", E_H_field, " (H-Field formulation)"
      else
         write(ui,*) "MAIN (B): action indef. avec E_H_field = ", E_H_field
         write(ui,*) "Aborting..."
         errco = -20
         return
      endif
      write(26,*) "   bloch_vec = ", bloch_vec
      write(26,*) "bloch_vec/pi = ", (bloch_vec(i)/D_PI,i=1,2)
      z_tmp = sqrt(shift_ksqr)/(2.0d0*D_PI)
      write(26,*) "shift_ksqr = ", shift_ksqr, z_tmp
      write(26,*) "integer super-vector :"
      write(26,*) "int_used, int_max, int_used/int_max   = ", int_used , int_max, dble(int_used)/dble(int_max)
      write(26,*) "cmplx super-vector : "
      write(26,*) "cmplx_used, cmplx_max, cmplx_used/cmplx_max = ", cmplx_used, cmplx_max, dble(cmplx_used)/dble(cmplx_max)
      write(26,*) "Real super-vector : "
      write(26,*) "real_used, real_max, real_max/real_used = ", real_used, real_max, dble(real_max)/dble(real_used)
      write(26,*)
      write(26,*) "n_modes, nvect, n_conv = ", n_modes, nvect, n_conv
      write(26,*) "nonz, n_msh_pts*n_modes, ", "nonz/(n_msh_pts*n_modes) = ", nonz, &
         n_msh_pts*n_modes, dble(nonz)/dble(n_msh_pts*n_modes)
      write(26,*) "nonz, nonz_max, nonz_max/nonz = ", nonz, nonz_max, dble(nonz_max)/dble(nonz)
      write(26,*) "nonz, int_used, int_used/nonz = ", nonz, int_used, dble(int_used)/dble(nonz)
!
!  write(26,*) "len_skyl, n_msh_pts*n_modes, len_skyl/(n_msh_pts*n_modes) = ",
!  *   len_skyl, n_msh_pts*n_modes, dble(len_skyl)/dble(n_msh_pts*n_modes)
!
      write(26,*)
      do i=1,n_modes
         write(26,"(i4,2(g22.14),g18.10)") i, v_eigs_beta(i)
      enddo
      write(26,*)
      write(26,*) "n_core = ", n_core
      write(26,*) "eps_eff = ", (eps_eff(i),i=1,n_typ_el)
      write(26,*) "v_refindex_n = ", (v_refindex_n(i),i=1,n_typ_el)
      write(26,*)
      write(26,*) "conjugate pair problem", pair_warning, "times"
      write(26,*)
      write(26,*) "mesh_file = ", mesh_file
      write(26,*) "gmsh_file = ", gmsh_file
      write(26,*) "log_file  = ", log_file
      close(26)
!
      write(ui,*) "   .   .   ."
      write(ui,*) "   .   .   ."
      write(ui,*) "   .   .   ."
      write(ui,*) "  and   we're  done!"
   endif

   write(ui,*) "-----------------------------------------------"
   write(ui,*)
!
   deallocate(a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L)

end subroutine calc_EM_modes

