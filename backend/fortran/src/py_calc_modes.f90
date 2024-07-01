#include "numbat_decl.h"

!  Solves the electromagnetic FEM problem defined in
!  Dossou & Fontaine, Comp Meth. App. Mech. Eng, 194, 837 (2005).

!  The weak formulation of Maxwell wave equation is in Eqs 14, 15.
!  \langle 1/\mu (\nabla_t \times E_t), (\nabla_t \times F_t) \rangle
!  - \omega^2 \langle (\epsilon E_t, F_t)
!  = \beta^2 \langle 1/\mu (\nabla_t hE_z -E_t, F_t), \rangle

!  \langle 1/\mu E_t, \nabla_t F_z \rangle
!  - \langle 1/\mu\nabla_t hE_z, \nabla_t F_z \rangle
!  + \omega^2 \langle\eps hE_z, F_z\rangle =0

!  where \hE_z = -1/\beta E_z

!  The fields are expanded in in-plane vector and longitudinal scalar elements
!  \vecphi_h and \psi_h:
!  E = E_{t,h} \vecphi_h + \unitz hE_{z,h} \psi_h = [E_{t,h} \vecphi_h, hE_{z,h} \psi_h ]
!  F = F_{t,h} \vecphi_h + \unitz F_{z,h} \psi_h   (note F, not hF)

!  Then  inner product (L_1 E, L_2 F) is evaluted:
!  (E,F) = \int dx dy   (L_2 F)^* \cdot (L_1 E)
!  = \int dx dy   ((L_2 F)_t)^* \cdot ((L_1 E)_t)
!  +  ((L_2 F)_z)^* . ((L_1 E)_z)

!  = \int dx dy   ((L_2 F)_t)^* \cdot ((L_1 E)_t)
!  +  ((L_2 F)_z)^* . ((L_1 E)_z)

!  This translates to the geneig problem (eq 40)

!  [ K_tt   0 ] [ E_t,h]  = \beta^2  [M_tt   (K_zt)^T] [E_t,h]
!  [ 0      0 ] [ hE_z,h]            [K_zt    K_zz   ] [hE_z,h]




!  lambda - free space wavelength in m
!  n_modes - desired number of eigenvectors
!  n_msh_pts - number of FEM mesh points
!  n_msh_el  - number of FEM (triang) elements
!  n_typ_el  - number of types of elements (and therefore elements)
!  v_refindex_n - array of effective index of materials
!  bloch_vec - in-plane k-vector (normally tiny just to avoid degeneracies)
!  shift_ksqr   - k_est^2 = n^2 vacwavenum_k0^2  : estimate of eigenvalue k^2
!  bnd_cnd_i - bnd conditions (Dirichlet = 0, Neumann = 1, Periodic = 2)
!  v_evals_beta_adj  - array of eigenvalues kz
!  m_evecs_adj   - 4-dim array of solutions [field comp, node of element (1..13)?!, eigvalue, element number] (strange ordering)
!  mode_pol  - unknown - never used in python
!  table_nod - 2D array [node_on_elt-1..6][n_msh_el] giving the mesh point mp of each node
!  Points where type_el[mp] is not the same for all 6 nodes must be interface points
!  type_el   - n_msh_el array: material index for each element
!  type_nod  - is boundary node?
!  mesh_xy  - (2 , n_msh_pts)  x,y coords?
!  ls_material  - (1, nodes_per_el+7, n_msh_el)

module calc_em_impl

   use numbatmod
   use class_stopwatch

contains

   subroutine calc_em_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &
      E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_typ_el, v_refindex_n, &
      v_evals_beta_adj, m_evecs_adj, mode_pol, table_nod, type_el, type_nod, mesh_xy, ls_material, errco, emsg)

      implicit none

      integer(8), intent(in) :: n_modes
      double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
      complex(8), intent(in) :: shift_ksqr

      integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
      character(len=*), intent(in) :: mesh_file
      integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_typ_el

      complex(8), intent(in) ::  v_refindex_n(n_typ_el)

      complex(8), target, intent(out) :: v_evals_beta_adj(n_modes)
      complex(8), target, intent(out) :: m_evecs_adj(3,nodes_per_el+7,n_modes,n_msh_el)

      complex(8), intent(out) :: mode_pol(4,n_modes)
      integer(8), intent(out) :: table_nod(nodes_per_el, n_msh_el)
      integer(8), intent(out) :: type_el(n_msh_el), type_nod(n_msh_pts)
      double precision, intent(out) :: mesh_xy(2,n_msh_pts)
      complex(8), intent(out) :: ls_material(1,nodes_per_el+7,n_msh_el)

      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg


      integer(8) neq


      integer(8) int_max, cmplx_max, int_used, cmplx_used
      integer(8) real_max, real_used

      integer(8), dimension(:), allocatable :: a_iwork
      complex(8), dimension(:), allocatable :: b_zwork
      double precision, dimension(:), allocatable :: c_dwork
      double precision, dimension(:,:), allocatable :: d_dwork
      double precision, dimension(:), allocatable :: e_dwork  !  take over work from b_zwork but have same shape

      integer(8), dimension(:), allocatable :: iindex
      complex(8), dimension(:,:), allocatable :: overlap_L

!  Declare the pointers of the integer super-vector
      integer(8) ip_table_E, ip_table_N_E_F, ip_visited
      integer(8) ip_type_N_E_F, ip_eq
      integer(8) ip_period_N, ip_nperiod_N
      integer(8) ip_period_N_E_F, ip_nperiod_N_E_F

!  Declare the pointers of the real super-vector
      integer(8) jp_x_N_E_F


      integer(8) jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer(8) jp_trav, jp_evecs
      complex(8) pp(n_typ_el), qq(n_typ_el)
      complex(8) eps_eff(n_typ_el)

      integer(8) n_msh_pts_p3, ui_out

!  Variable used by valpr
      integer(8) dim_krylov, ltrav
      integer(8) n_conv, i_base
      !double precision ls_data(10)

      integer(8) n_core(2)  !  index of highest epsilon material, seems funky
      integer(8) n_edge, n_face, n_ddl, n_ddl_max


      double precision vacwavenum_k0, dim_x, dim_y

      double precision time_fact, time_arpack

!  Declare the pointers of the real super-vector
      integer(8) kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im
      integer(8) kp_mat1_re, kp_mat1_im

!  Declare the pointers of for sparse matrix storage
      integer(8) ip_col_ptr, ip_row
      integer(8) jp_mat2
      integer(8) ip_work, ip_work_sort, ip_work_sort2
      integer(8) nonz, nonz_max, max_row_len


      complex(8), pointer :: p_sol (:,:,:,:)
      complex(8), pointer :: p_beta(:)


      integer(8) :: ilo, ihi, i_md

      integer :: is_em, alloc_stat, alloc_remote
      double precision tol

      type(Stopwatch) :: clock_main, clock_spare


      is_em = 1
      alloc_remote = 0
      ui_out = stdout


      !  TODO: unallocated arrays can not be passed as function arguments
      !  Can be done by quasi-globals variables in a module

      if (alloc_remote > 0) then
         !call prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, int_max, cmplx_max, real_max, &
         !  a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L,  errco, emsg)
         !RETONERROR(errco)
         write(ui_out,*) 'alloc remote in calc_em_impl is broken'
         call exit(1)
      else

         call array_size(n_msh_pts, n_msh_el, n_modes, &
            int_max, cmplx_max, real_max, n_ddl, errco, emsg)
         RETONERROR(errco)

         allocate(a_iwork(int_max), STAT=alloc_stat)
         call check_alloc(alloc_stat, int_max, "a_iwork", 101, errco, emsg)
         RETONERROR(errco)

         allocate(b_zwork(cmplx_max), STAT=alloc_stat)
         call check_alloc(alloc_stat, cmplx_max, "b_zwork", 101, errco, emsg)
         RETONERROR(errco)

         allocate(c_dwork(real_max), STAT=alloc_stat)
         call check_alloc(alloc_stat, real_max, "c_dwork", 101, errco, emsg)
         RETONERROR(errco)

         allocate(iindex(n_modes), STAT=alloc_stat)
         call check_alloc(alloc_stat, n_modes, "iindex", 101, errco, emsg)
         RETONERROR(errco)

         if (is_em > 0) then
            allocate(d_dwork(2,n_ddl), STAT=alloc_stat)
            call check_alloc(alloc_stat, 2*n_ddl, "d_dwork", 101, errco, emsg)
            RETONERROR(errco)

            allocate(e_dwork(cmplx_max), STAT=alloc_stat)
            call check_alloc(alloc_stat, cmplx_max, "e_dwork", 101, errco, emsg)
            RETONERROR(errco)

            allocate(overlap_L(n_modes,n_modes), STAT=alloc_stat)
            call check_alloc(alloc_stat, n_modes*n_modes, "overlap_L", 101, errco, emsg)
            RETONERROR(errco)

         endif
      endif



!  nsym = 1 !  nsym = 0 => symmetric or hermitian matrices


      dim_krylov = 2*n_modes + n_modes/2 +3

      call clock_main%reset()


      dim_x = dimscale_in_m
      dim_y = dimscale_in_m

      !  Fill:  mesh_xy, type_nod, type_el, table_nod
      call construct_fem_node_tables (n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, dim_x, dim_y, mesh_file, &
         mesh_xy, type_nod, type_el, table_nod, errco, emsg)
      RETONERROR(errco)

      !  Storage locations in sequence
      !  - table_edge_face = a_iwork(ip_table_N_E_F),   shape: 14 x n_msh_el
      !  - visited         = a_iwork(ip_visited),       shape: n_ddl_max = npt + n_msh_el = 4 n_msh_el
      !  - table_edges     = a_iwork(ip_table_E)        shape: 4 x n_msh_pts
      !
      !  visited is used as workspace. has no meaning between functions
      !
      !  V = number of vertices
      !  E = number of edges
      !  F = number of faces
      !  C = number of cells (3D, tetrahedron)
      !
      !  From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
      !  n_msh_pts = (number of vertices) + (number of mid-edge point) = V + E;
      !
      !  neq and nonz are some kind of dimension for the left and right eigenoperators

      !  TODO: move next three calls into a single  construct_table_N_E_F procedure

      !  Fills:  table_edge_face[1,:]
      ip_table_N_E_F = 1
      call list_face (n_msh_el, a_iwork(ip_table_N_E_F))

      !  n_ddl_max = max(N_Vertices) + max(N_Edge) + max(N_Face)
      !  For P2 FEM n_msh_pts=N_Vertices+N_Edge
      !  note: each element has 1 face, 3 edges and 10 P3 nodes
      !  so table_N_E_F = table_edge_face has dimensions 14 x n_msh_el

      !  each element is a face
      n_face = n_msh_el

      n_ddl_max = n_msh_pts + n_face

      ip_visited =  ip_table_N_E_F  + 14*n_msh_el
      ip_table_E = ip_visited + n_ddl_max

      !  Fills: n_edge, table_edge[1..4,:], table_edge_face[2:4,:], visited[1:n_msh_pts]
      !  Todo!  move n_edge later in list as an out variable
      call list_edge (n_msh_el, n_msh_pts, nodes_per_el, n_edge, type_nod, table_nod, &
         a_iwork(ip_table_E), a_iwork(ip_table_N_E_F), a_iwork(ip_visited))

      !  Fills: remainder of table_edge_face[5:,:], visited[1:n_msh_pts], n_msh_pts_3
      !  Todo: move n_msh_pts_p3 later
      call list_node_P3 (n_msh_el, n_msh_pts, nodes_per_el, n_edge, n_msh_pts_p3, table_nod, &
         a_iwork(ip_table_N_E_F), a_iwork(ip_visited))

      !  TODO: what is signif of this quanitty?
      n_ddl = n_edge + n_face + n_msh_pts_p3


      if (debug .eq. 1) then
         write(ui_out,*) "py_calc_modes.f: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
         write(ui_out,*) "py_calc_modes.f: n_msh_pts_p3 = ", n_msh_pts_p3
         write(ui_out,*) "py_calc_modes.f: n_vertex, n_edge, n_face,", " n_msh_el = ", &
            (n_msh_pts - n_edge), n_edge, n_face, n_msh_el
         write(ui_out,*) "py_calc_modes.f: 2D case of the Euler &
         & characteristic: V-E+F=1-(number of holes)"
         write(ui_out,*) "py_calc_modes.f: Euler characteristic: V - E + F &
         &= ", (n_msh_pts - n_edge) - n_edge + n_face
      endif


!C  overwriting pointers ip_row_ptr, ..., ip_adjncy

      ip_type_N_E_F = ip_table_E + 4*n_edge   !  not sure why 4* n_edge, not 4*n_msh_pts?


      !  TODO:
      !  ip is an index into an a_iwork, make this clearer!
      !  jp is an index into an b_zwork
      !  kp is an index into an c_dwork

      !  Offsets into the b_zwork workspace
      jp_x_N_E_F = 1


      !  Offsets into the a_iwork workspace
      ip_period_N = ip_type_N_E_F + 2*n_ddl
      ip_nperiod_N = ip_period_N + n_msh_pts
      ip_period_N_E_F = ip_nperiod_N + n_msh_pts
      ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
      ip_eq = ip_nperiod_N_E_F + n_ddl


      !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
      !  Should be using c_dwork for x_E_F ?
      call type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, type_nod, table_nod, &
         a_iwork(ip_table_N_E_F), a_iwork(ip_visited), a_iwork(ip_type_N_E_F), mesh_xy, &
      !b_zwork(jp_x_N_E_F) &
         d_dwork &
         )


      !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
      call get_coord_p3 (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, table_nod, type_nod, &
         a_iwork(ip_table_N_E_F), a_iwork(ip_type_N_E_F), mesh_xy, &
      !b_zwork(jp_x_N_E_F), &
         d_dwork, &
         a_iwork(ip_visited))




      !  TODO: the b_zwork should actually be the d_dwork containing x_N_E_F, but only matters for periodic
      call set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el, mesh_xy, nodes_per_el, &
         type_nod, table_nod, n_ddl, neq, ip_type_N_E_F, ip_eq, a_iwork, &
      !b_zwork, &
         d_dwork,  &
         int_max, debug)


!  Needed vars from above here:  ip_eq, jp_x_N_E_F, ip_period_N


!  Sparse matrix storage

      ip_col_ptr = ip_eq + 3*n_ddl

      call csr_max_length (n_msh_el, n_ddl, neq, a_iwork(ip_table_N_E_F), &
         a_iwork(ip_eq), a_iwork(ip_col_ptr), nonz_max)

!  ip = ip_col_ptr + neq + 1 + nonz_max
         !ip = ip_col_ptr + neq + 1
         ip_row = ip_col_ptr + neq + 1

      if (ip_row .gt. int_max) then
         write(emsg,*) "py_calc_modes.f: ip_row > int_max : ", ip_row, int_max, "py_calc_modes.f: nonz_max = ", &
            nonz_max, "py_calc_modes.f: increase the size of int_max"
         errco = -11
         return
      endif


      call csr_length (n_msh_el, n_ddl, neq,  a_iwork(ip_table_N_E_F), a_iwork(ip_eq), a_iwork(ip_row), &
         a_iwork(ip_col_ptr), nonz_max, nonz, max_row_len, ip_row, int_max, debug)

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*n_ddl
      ip_work_sort2 = ip_work_sort + max_row_len

!  sorting csr ...
      call sort_csr (neq, nonz, max_row_len, a_iwork(ip_row), a_iwork(ip_col_ptr), a_iwork(ip_work_sort), a_iwork(ip_work), &
         a_iwork(ip_work_sort2))

      if (debug .eq. 1) then
         write(ui_out,*) "py_calc_modes.f: nonz_max = ", nonz_max
         write(ui_out,*) "py_calc_modes.f: nonz = ", nonz
         write(ui_out,*) "py_calc_modes.f: cmplx_max/nonz = ", dble(cmplx_max)/dble(nonz)
      endif

      int_used = ip_work_sort2 + max_row_len

      if (int_max .lt. int_used) then
         write(emsg,*)'The size of the integer supervector is too small', 'integer super-vec: int_max  = ', &
            int_max, 'integer super-vec: int_used = ', int_used
         errco = -12
         return
      endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      jp_mat2 = jp_x_N_E_F + 3*n_ddl

      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq

!  Eigenvectors
      jp_vschur = jp_resid + neq
      jp_trav = jp_vschur + neq*dim_krylov

      ltrav = 3*dim_krylov*(dim_krylov+2)
      jp_evecs = jp_trav + ltrav

      cmplx_used = jp_evecs + neq*n_modes

      if (cmplx_max .lt. cmplx_used)  then
         write(emsg,*)'The size of the complex supervector is too small', 'complex super-vec: int_max  = ', &
            cmplx_max, 'complex super-vec: int_used = ', cmplx_used
         errco = -13
         return
      endif

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


!###############################################

!  ----------------------------------------------------------------
!  convert from 1-based to 0-based
!  ----------------------------------------------------------------

!  do j = 1, neq+1
!  a_iwork(j+ip_col_ptr-1) = a_iwork(j+ip_col_ptr-1) - 1
!  end do
!  do  j = 1, nonz
!  a_iwork(j+ip_row-1) = a_iwork(j+ip_row-1) - 1
!  end do
!

      ilo = ip_col_ptr-1 + 1
      ihi = ip_col_ptr-1 + neq + 1
      a_iwork(ilo:ihi) = a_iwork(ilo:ihi) - 1

      ilo = ip_row-1 + 1
      ihi = ip_row-1 + nonz
      a_iwork(ilo:ihi) = a_iwork(ilo:ihi) - 1




!  The CSC indexing, i.e., ip_col_ptr, is 1-based
!  (but valpr.f will change the CSC indexing to 0-based indexing)
      i_base = 0


      write(ui_out,*)
      write(ui_out,*) "-----------------------------------------------"


      vacwavenum_k0 = 2.0d0*D_PI/lambda


      call  check_materials_and_fem_formulation(E_H_field,n_typ_el, &
         vacwavenum_k0, v_refindex_n, eps_eff, n_core, pp, qq, debug, ui_out, errco, emsg)
      RETONERROR(errco)


!  Main eigensolver
      write(ui_out,*) "EM FEM: "

      p_sol  => m_evecs_adj
      p_beta => v_evals_beta_adj


!  Assemble the coefficient matrix A and the right-hand side F of the
!  finite element equations

      write(ui_out,'(A,A)') "   - assembling linear system "

      call clock_spare%reset()


      !  Build the actual matrices A (mat_1) and M(mat_2) for the arpack solving.  (M = identity?)
      call asmbly (bdy_cdn, i_base, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, &
         shift_ksqr, bloch_vec, n_typ_el, pp, qq, &
         table_nod, a_iwork(ip_table_N_E_F), type_el, &
         a_iwork(ip_eq), a_iwork(ip_period_N), a_iwork(ip_period_N_E_F), &
         mesh_xy, &
      !b_zwork(jp_x_N_E_F), &
         d_dwork, &
         nonz, a_iwork(ip_row), a_iwork(ip_col_ptr), &
         c_dwork(kp_mat1_re), c_dwork(kp_mat1_im), b_zwork(jp_mat2), a_iwork(ip_work))

      call clock_spare%stop()
      write(ui_out,'(A,A)') '      ', clock_spare%to_string()

      write(ui_out,*) "  - solving linear system"

      call clock_spare%reset()

      !  This is the main solver.
      !  On completion:
      !  unshifted unsorted eigenvalues are in p_beta[1..n_modes]
      !  eigvectors are in are b_zwork[jp_evecs..?]

      !  TODO: following are no longer needed:  b_zwork(jp_trav/vect1/vect2),

      call valpr_64 (i_base, &
      !b_zwork(jp_vect1), &  !  unused
      !b_zwork(jp_vect2), &  !  unused
      !b_zwork(jp_trav), &  !  unused
         !b_zwork(jp_workd), b_zwork(jp_resid), ltrav,  &
         dim_krylov, n_modes, neq, itermax,  &
         tol, nonz, a_iwork(ip_row), a_iwork(ip_col_ptr), &
         c_dwork(kp_mat1_re), c_dwork(kp_mat1_im), b_zwork(jp_mat2), &
         c_dwork(kp_lhs_re), c_dwork(kp_lhs_im), c_dwork(kp_rhs_re), c_dwork(kp_rhs_im), &
         p_beta, b_zwork(jp_vschur), b_zwork(jp_evecs), &
         n_conv, time_fact, time_arpack, debug, errco, emsg)
      RETONERROR(errco)



      if (n_conv .ne. n_modes) then
         write(emsg,*) "Convergence problem in valpr_64: n_conv != n_modes : ", &
            n_conv, n_modes ,"You should probably increase resolution of mesh!"
         errco = -19
         return
      endif

      call clock_spare%stop()
      write(ui_out,'(A,A)') '      ', clock_spare%to_string()


      call rescale_and_sort_eigensolutions(n_modes, shift_ksqr, p_beta, iindex)


      call clock_spare%reset()

      write(ui_out,*) "  - assembling eigen solutions"



      !  The eigenvectors will be stored in the array sol
      !  The eigenvalues and eigenvectors are renumbered
      !  using the permutation vector iindex
      call array_sol ( bdy_cdn, n_modes, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, &
         n_core, bloch_vec, iindex, table_nod, a_iwork(ip_table_N_E_F), type_el, &
         a_iwork(ip_eq), a_iwork(ip_period_N), a_iwork(ip_period_N_E_F), &
         mesh_xy, &
      !b_zwork(jp_x_N_E_F),
         !d_dwork, &  !  this should be an e_ework
         e_dwork, &  !  this should be an e_ework
         p_beta, mode_pol, b_zwork(jp_evecs), p_sol , errco, emsg)
      RETONERROR(errco)


      write(ui_out,*) "  - finding mode energies "
      !  Calculate energy in each medium (typ_el)
      call mode_energy (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_core, &
         table_nod, type_el, n_typ_el, eps_eff,&
         mesh_xy, p_sol, p_beta, mode_pol)




      !  Doubtful that this check is of any value: delete?
      !  call check_orthogonality_of_em_sol(n_modes, n_msh_el, n_msh_pts, n_typ_el, pp, table_nod, &
      !  type_el, mesh_xy, v_evals_beta_adj, m_evecs_adj, &!v_evals_beta_pri, m_evecs_pri,
      !  overlap_L, overlap_file, debug, ui_out, &
      !  pair_warning, vacwavenum_k0, errco, emsg)
      !  RETONERROR(errco)


      !  Should this happen _before_ check_ortho?


      !  The z-component must be multiplied by -ii*beta in order to
      !  get the physical, un-normalised z-component
      !  (see Eq. (25) of the JOSAA 2012 paper)
      !  TODO: is this really supposed to be x i beta , or just x beta  ?
      do i_md=1,n_modes
         !!do iel=1,n_msh_el
         !!  m_evecs_adj(3,inod,i_md,iel) = C_IM_ONE * p_beta(i_md) * m_evecs_adj(3,inod,i_md,iel)
         !!enddo

         !do inod=1,nodes_per_el+7
         !  m_evecs_adj(3,inod,i_md,:) = C_IM_ONE * p_beta(i_md) * m_evecs_adj(3,inod,i_md,:)
         !enddo

         m_evecs_adj(3,:,i_md,:) = C_IM_ONE * p_beta(i_md) * m_evecs_adj(3,:,i_md,:)

      enddo


      call array_material_EM (n_msh_el, n_typ_el, v_refindex_n, type_el, ls_material)

      !  Normalisation. Can't use this if we don't do check_ortho.  Not needed
      !  call normalise_fields(n_modes, n_msh_el, nodes_per_el, m_evecs_adj, overlap_L)

      write(ui_out,*) "  - finished"
      !if (debug .eq. 1) then
      !  write(ui_out,*) "py_calc_modes.f: CPU time for normalisation :", (time2_J-time1_J)
      !endif
      !
      !  Orthonormal integral
      !  if (debug .eq. 1) then
      !  write(ui_out,*) "py_calc_modes.f: Product of normalised field"
      !  overlap_file = "Orthogonal_n.txt"
      !  call get_clocks( systime1_J, time1_J)
      !  call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, pp, table_nod, &
      !  type_el, mesh_xy, v_evals_beta_adj, v_evals_beta_pri, m_evecs_adj, m_evecs_pri, overlap_L, overlap_file, debug, &
      !  pair_warning, vacwavenum_k0)
      !  call get_clocks( systime2_J, time2_J)
      !  write(ui_out,*) "py_calc_modes.f: CPU time for orthogonal :", (time2_J-time1_J)
      !  endif
      !


      call clock_spare%stop()
      call clock_main%stop()

      deallocate(a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L)

      write(ui_out,*) "-----------------------------------------------"

      !  call report_results_em(debug, ui_out, &
      !  n_msh_pts, n_msh_el, &
      !  time1, time2, time_fact, time_arpack,  time1_postp, time2_postp, &
      !  lambda, e_h_field, bloch_vec, bdy_cdn,  &
      !  int_max, cmplx_max, cmplx_used,  n_core, n_conv, n_modes, &
      !  n_typ_el, neq, nonz_max, dim_krylov, &
      !  shift_ksqr, v_evals_beta_adj, eps_eff, v_refindex_n)



   end subroutine calc_em_modes_impl

   subroutine check_materials_and_fem_formulation(E_H_field,n_typ_el, &
      vacwavenum_k0, v_refindex_n, eps_eff, n_core, pp, qq, debug, ui_out, errco, emsg)

      integer(8), intent(in) :: E_H_field, debug
      integer(8), intent(in) :: n_typ_el, ui_out
      double precision, intent(in):: vacwavenum_k0
      complex(8), intent(in) :: v_refindex_n(n_typ_el)
      complex(8), intent(out) :: eps_eff(n_typ_el)

      integer(8), intent(out) :: n_core(2)
      complex(8), intent(out) :: pp(n_typ_el), qq(n_typ_el)
      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) i
      logical is_homogeneous

      eps_eff = v_refindex_n**2

      !  what actually even is this?
      if(dble(eps_eff(1)) .gt. dble(eps_eff(2))) then
         n_core(1) = 1
      else
         n_core(1) = 2
      endif
      n_core(2) = n_core(1)

!  Check that the structure is not entirely homogeneous (TODO: does this actually matter?)
      is_homogeneous = .true.
      do i=1,n_typ_el-1

         if (.not. almost_equal(dble(eps_eff(i)), dble(eps_eff(i+1)))) then
            is_homogeneous = .false.
         elseif (.not. almost_equal(dimag(eps_eff(i)), dimag(eps_eff(i+1)))) then
            is_homogeneous = .false.
         endif

      enddo

      if (is_homogeneous) then
         emsg = "py_calc_modes.f: FEM routine cannot adjacent identical layers. Define layer as object.ThinFilm."
         errco = -17
         return
      endif


      if(debug .eq. 1) then
         write(ui_out,*) "py_calc_modes.f: n_core = ", n_core
         if(E_H_field .eq. FEM_FORMULATION_E) then
            write(ui_out,*) "py_calc_modes.f: E-Field formulation"
         else
            write(ui_out,*) "py_calc_modes.f: H-Field formulation"
         endif
      endif


      !  set up some kind of mass vectors for the FEM
      !  weird place but ok.
      if(E_H_field .eq. FEM_FORMULATION_E) then
         qq = eps_eff*vacwavenum_k0**2
         pp = 1.0d0
      elseif(E_H_field .eq. FEM_FORMULATION_H) then
         qq = vacwavenum_k0**2
         pp = 1.0d0/eps_eff
      endif

   end subroutine

   subroutine check_orthogonality_of_em_sol(n_modes, n_msh_el, n_msh_pts, n_typ_el, pp, table_nod, &
      type_el, mesh_xy, v_evals_beta_adj, m_evecs_adj, &
   !v_evals_beta_pri, m_evecs_pri, &
      overlap_L, overlap_file, debug, ui_out, pair_warning, vacwavenum_k0, errco, emsg)

      use numbatmod
      logical pair_warning

      integer(8), intent(in) :: n_modes, debug, ui_out
      integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_typ_el
      complex(8) pp(n_typ_el)

      integer(8), intent(out) :: table_nod(nodes_per_el, n_msh_el)
      integer(8), intent(out) :: type_el(n_msh_el)
      double precision, intent(out) :: mesh_xy(2,n_msh_pts)
      double precision vacwavenum_k0

      complex(8), target, intent(out) :: v_evals_beta_adj(n_modes)
      complex(8), target, intent(out) :: m_evecs_adj(3,nodes_per_el+7,n_modes,n_msh_el)

      complex(8), dimension(:,:) :: overlap_L


      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      character(len=FNAME_LENGTH)  overlap_file

      !complex(8)  :: v_evals_beta_pri(n_modes)
      !complex(8)  :: m_evecs_pri(3,nodes_per_el+7,n_modes,n_msh_el)

      !  Orthogonal integral
      pair_warning = .false.

      if (debug .eq. 1) then
         write(ui_out,*) "py_calc_modes.f: Field product"
      endif

      overlap_file = "Orthogonal.txt"

      call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_typ_el, pp, table_nod, &
         type_el, mesh_xy, v_evals_beta_adj, m_evecs_adj, &
      !v_evals_beta_pri, m_evecs_pri,
         overlap_L, overlap_file, debug, pair_warning, vacwavenum_k0)

      if (pair_warning .and. n_modes .le. 20) then
         emsg = "py_calc_modes.f: Warning found 1 BM of cmplx conj pair, increase num_BMs to include the other."
         errco = -57
      endif


   end subroutine

   subroutine report_results_em(debug, ui_out, &
      n_msh_pts, n_msh_el, &
      time1, time2, time_fact, time_arpack, time1_postp, &
      lambda, e_h_field, bloch_vec, bdy_cdn,  &
      int_max, cmplx_max, cmplx_used,  n_core, n_conv, n_modes, &
      n_typ_el, neq, nonz_max, dim_krylov, &
      shift_ksqr, v_evals_beta_adj, eps_eff, v_refindex_n)



      use numbatmod

      integer(8) debug, ui_out, e_h_field, bdy_cdn
      integer(8) int_max, cmplx_max, cmplx_used, int_used, real_max, real_used, n_msh_pts, n_msh_el
      double precision bloch_vec(2), lambda
      double precision time1, time2, start_time, end_time, time_fact, time_arpack, time1_postp
      integer(8) n_conv, n_modes, n_typ_el, nonz, nonz_max, n_core(2), neq, dim_krylov
      character(len=FNAME_LENGTH)  log_file
      complex(8), intent(in) :: shift_ksqr
      complex(8), target, intent(out) :: v_evals_beta_adj(n_modes)
      complex(8) eps_eff(n_typ_el)

      complex(8), intent(in) ::  v_refindex_n(n_typ_el)

      complex(8) z_tmp
      integer(8) i

      !  TODO: hook these up if needed
      cmplx_max = 0
      int_max = 0
      int_used =0
      nonz = 0
      nonz_max = 0
      cmplx_used = 0
      real_max = 0
      real_used = 0

      if (debug .eq. 1) then
         write(ui_out,*)
         write(ui_out,*) 'Total CPU time (sec.)  = ', (time2-time1)
         open (unit=26,file=log_file)
         write(26,*)
         write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
         write(26,*) "Start time   = ", start_time
         write(26,*) "End time  = ", end_time
         write(26,*) "Total CPU time (sec.) = ", (time2-time1)
         write(26,*) "LU factorisation : CPU time and % Total time = ", time_fact, &
            100*(time_fact)/(time2-time1),"%"
         write(26,*) "ARPACK : CPU time and % Total time = ", time_arpack, &
            100*(time_arpack)/(time2-time1),"%"
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
         endif
         write(26,*) "   bloch_vec = ", bloch_vec
         write(26,*) "bloch_vec/pi = ", (bloch_vec(i)/D_PI,i=1,2)
         z_tmp = sqrt(shift_ksqr)/(2.0d0*D_PI)
         write(26,*) "shift_ksqr = ", shift_ksqr, z_tmp
         !  write(26,*) "integer super-vector :"
         !  write(26,*) "int_used, int_max, int_used/int_max   = ", int_used , int_max, dble(int_used)/dble(int_max)
         !write(26,*) "cmplx super-vector : "
         !write(26,*) "cmplx_used, cmplx_max, cmplx_used/cmplx_max = ", cmplx_used, cmplx_max, dble(cmplx_used)/dble(cmplx_max)

         !write(26,*) "Real super-vector : "
         !write(26,*) "real_used, real_max, real_max/real_used = ", real_used, real_max, dble(real_max)/dble(real_used)
         write(26,*)
         write(26,*) "n_modes, dim_krylov, n_conv = ", n_modes, dim_krylov, n_conv
         !write(26,*) "nonz, n_msh_pts*n_modes, ", "nonz/(n_msh_pts*n_modes) = ", nonz, &
         !  n_msh_pts*n_modes, dble(nonz)/dble(n_msh_pts*n_modes)
         !write(26,*) "nonz, nonz_max, nonz_max/nonz = ", nonz, nonz_max, dble(nonz_max)/dble(nonz)
         !write(26,*) "nonz, int_used, int_used/nonz = ", nonz, int_used, dble(int_used)/dble(nonz)

!  write(26,*) "len_skyl, n_msh_pts*n_modes, len_skyl/(n_msh_pts*n_modes) = ",
!  *   len_skyl, n_msh_pts*n_modes, dble(len_skyl)/dble(n_msh_pts*n_modes)

         write(26,*)
         do i=1,n_modes
            write(26,"(i4,2(g22.14),g18.10)") i, v_evals_beta_adj(i)
         enddo
         write(26,*)
         write(26,*) "n_core = ", n_core
         write(26,*) "eps_eff = ", (eps_eff(i),i=1,n_typ_el)
         write(26,*) "v_refindex_n = ", (v_refindex_n(i),i=1,n_typ_el)
         write(26,*)
         !write(26,*) "conjugate pair problem", pair_warning, "times"
         write(26,*)
         !write(26,*) "mesh_file = ", mesh_file
         !write(26,*) "gmsh_file = ", gmsh_file
         !write(26,*) "log_file  = ", log_file
         close(26)

      endif

      write(ui_out,*) "-----------------------------------------------"
      write(ui_out,*)

   end subroutine




   subroutine rescale_and_sort_eigensolutions(n_modes, shift_ksqr, p_beta, iindex)

      integer(8), intent(in) :: n_modes
      complex(8), intent(in) :: shift_ksqr
      complex(8), pointer :: p_beta(:)
      integer(8), dimension(:), allocatable :: iindex

      integer(8) i

      complex(8) z_beta

      !TODO: make a function. Turn beta^2 raw eig into actual beta
      do i=1,n_modes
         !  z_tmp0 = p_beta(i)
         !  z_tmp = 1.0d0/z_tmp0+shift_ksqr
         !  z_beta = sqrt(z_tmp)

         z_beta = sqrt(1.0d0/p_beta(i)+shift_ksqr )

         !  Mode classification - we want the forward propagating mode
         if (abs(imag(z_beta)/z_beta) .lt. 1.0d-8) then
            !  re(z_beta) > 0 for forward propagating mode
            if (dble(z_beta) .lt. 0) z_beta = -z_beta
         else
            !  im(z_beta) > 0 for forward decaying evanescent mode  !rarely relevant for us
            if (imag(z_beta) .lt. 0) z_beta = -z_beta
         endif

         p_beta(i) = z_beta
      enddo



      !  order p_beta by magnitudes and store in iindex
      call z_indexx (n_modes, p_beta, iindex)
   end subroutine





end module calc_em_impl




