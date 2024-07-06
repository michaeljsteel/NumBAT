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
 !  n_elt_mats  - number of types of elements (and therefore elements)
 !  v_refindex_n - array of effective index of materials
 !  bloch_vec - in-plane k-vector (normally tiny just to avoid degeneracies)
 !  shift_ksqr   - k_est^2 = n^2 vacwavenum_k0^2  : estimate of eigenvalue k^2
 !  bnd_cnd_i - bnd conditions (Dirichlet = 0, Neumann = 1, Periodic = 2)
 !  v_evals_beta  - array of eigenvalues kz
 !  m_evecs   - 4-dim array of solutions [field comp, node of element (1..13)?!, eigvalue, element number] (strange ordering)
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
   use alloc

   use nbinterfaces
   use nbinterfacesb


contains

   subroutine calc_em_modes_impl( n_modes, lambda, dimscale_in_m, bloch_vec, shift_ksqr, &
      E_H_field, bdy_cdn, itermax, debug, mesh_file, n_msh_pts, n_msh_el, n_elt_mats, v_refindex_n, &
      v_evals_beta, m_evecs, mode_pol, table_nod, type_el, type_nod, mesh_xy, ls_material, errco, emsg)

      implicit none

      integer(8), parameter :: nodes_per_el = 6

      integer(8), intent(in) :: n_modes
      double precision, intent(in) :: lambda, dimscale_in_m, bloch_vec(2)
      complex(8), intent(in) :: shift_ksqr

      integer(8), intent(in) :: E_H_field, bdy_cdn, itermax, debug
      character(len=*), intent(in) :: mesh_file
      integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_elt_mats

      complex(8), intent(in) ::  v_refindex_n(n_elt_mats)

      complex(8), target, intent(out) :: v_evals_beta(n_modes)
      complex(8), target, intent(out) :: m_evecs(3,nodes_per_el+7,n_modes,n_msh_el)

      complex(8), intent(out) :: mode_pol(4,n_modes)
      integer(8), intent(out) :: table_nod(nodes_per_el, n_msh_el)
      integer(8), intent(out) :: type_el(n_msh_el), type_nod(n_msh_pts)
      double precision, intent(out) :: mesh_xy(2,n_msh_pts)
      complex(8), intent(out) :: ls_material(1,nodes_per_el+7,n_msh_el)

      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      !  ----------------------------------------------
      !  workspaces

      integer(8) neq,n_ddl, nonz

      !integer(8) n_edge, n_face, n_msh_pts_p3,


      integer(8) int_max, cmplx_max, real_max

      double precision, dimension(:,:), allocatable :: xy_N_E_F
      integer(8), dimension(:,:), allocatable :: table_N_E_F
      integer(8), dimension(:,:), allocatable :: type_N_E_F

      integer(8), dimension(:,:), allocatable :: m_eqs

      integer(8), dimension(:), allocatable :: v_eig_index
      complex(8), dimension(:,:), allocatable :: overlap_L

      complex(8), dimension(:,:), allocatable :: arp_evecs

      complex(8), dimension(:), allocatable :: mOp_stiff
      complex(8), dimension(:), allocatable :: mOp_mass

      integer(8), dimension(:), allocatable :: v_row_ind
      integer(8), dimension(:), allocatable :: v_col_ptr

      !integer(8), dimension(:), allocatable :: visited

      ! Currenly, periodic is not active

      integer(8), dimension(:), allocatable :: iperiod_N
      integer(8), dimension(:), allocatable :: iperiod_N_E_F
      integer(8), dimension(:), allocatable :: inperiod_N
      integer(8), dimension(:), allocatable :: inperiod_N_E_F


      ! Should these be dynamic?
      complex(8) pp(n_elt_mats), qq(n_elt_mats)
      complex(8) eps_eff(n_elt_mats)


      !  ----------------------------------------------


      integer(8) ui_out
      integer(8) :: i_md

      !  Variable used by valpr
      integer(8) dim_krylov
      integer(8) i_base
      double precision arp_tol



      integer(8) n_core(2)  !  index of highest epsilon material, seems funky
      double precision vacwavenum_k0, dim_x, dim_y



      type(Stopwatch) :: clock_main, clock_spare



      ui_out = stdout

      arp_tol = 1.0d-12 ! TODO: ARPACK_ stopping precision,  connect  to user switch

      call array_size(n_msh_pts, n_msh_el, n_modes, &
         int_max, cmplx_max, real_max, n_ddl, errco, emsg) !Only useful number out of here is n_ddl
      RETONERROR(errco)



      call double_alloc_2d(xy_N_E_F, 2_8, n_ddl, 'xy_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(type_N_E_F, 2_8, n_ddl, 'type_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(table_N_E_F, 14_8, n_msh_el, 'table_N_E_F', errco, emsg); RETONERROR(errco)

      call integer_alloc_2d(m_eqs, 3_8, n_ddl, 'm_eqs', errco, emsg); RETONERROR(errco)

      call integer_alloc_1d(v_eig_index, n_modes, 'v_eig_index', errco, emsg); RETONERROR(errco)
      call complex_alloc_2d(overlap_L, n_modes, n_modes, 'overlap_L', errco, emsg); RETONERROR(errco)

      call integer_alloc_1d(iperiod_N, n_msh_pts, 'iperiod_N', errco, emsg); RETONERROR(errco)
      call integer_alloc_1d(iperiod_N_E_F, n_ddl, 'iperiod_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_1d(inperiod_N, n_msh_pts, 'inperiod_N', errco, emsg); RETONERROR(errco)
      call integer_alloc_1d(inperiod_N_E_F, n_ddl, 'inperiod_N_E_F', errco, emsg); RETONERROR(errco)



      call clock_main%reset()


      dim_x = dimscale_in_m
      dim_y = dimscale_in_m

      !  Fill:  mesh_xy, type_nod, type_el, table_nod
      call construct_fem_node_tables (mesh_file, dim_x, dim_y, n_msh_el, n_msh_pts, nodes_per_el, n_elt_mats,   &
         mesh_xy, type_nod, type_el, table_nod, errco, emsg)
      RETONERROR(errco)


      call build_mesh_tables( &
      n_msh_el, n_msh_pts, nodes_per_el, n_ddl, &
      type_nod, table_nod, mesh_xy, &
      type_N_E_F, table_N_E_F, xy_N_E_F, &
      debug, errco, emsg)



      ! !  Storage locations in sequence
      ! !  - table_edge_face = table_N_E_F,   shape: 14 x n_msh_el
      ! !  - visited         shape: npt + n_msh_el = 4 n_msh_el
      ! !  - table_edges     shape: 4 x n_msh_pts
      ! !
      ! !  visited is used as workspace. has no meaning between functions
      ! !
      ! !  V = number of vertices
      ! !  E = number of edges
      ! !  F = number of faces
      ! !  C = number of cells (3D, tetrahedron)
      ! !
      ! !  From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
      ! !  n_msh_pts = (number of vertices) + (number of mid-edge point) = V + E;
      ! !
      ! !  neq and nonz are some kind of dimension for the left and right eigenoperators

      ! !  TODO: move next three calls into a single  construct_table_N_E_F procedure

      ! !  Fills:  table_edge_face[1,:]
      ! call list_face (n_msh_el, table_N_E_F)

      ! !  For P2 FEM n_msh_pts=N_Vertices+N_Edge
      ! !  note: each element has 1 face, 3 edges and 10 P3 nodes
      ! !  so table_N_E_F = table_edge_face has dimensions 14 x n_msh_el

      ! !  each element is a face
      ! n_face = n_msh_el

      ! !  Fills: n_edge, table_edge[1..4,:], table_edge_face[2:4,:], visited[1:n_msh_pts]
      ! !  Todo!  move n_edge later in list as an out variable
      ! call list_edge (n_msh_el, n_msh_pts, nodes_per_el, n_edge, type_nod, table_nod, &
      !    table_N_E_F, visited)

      ! !  Fills: remainder of table_edge_face[5:,:], visited[1:n_msh_pts], n_msh_pts_3
      ! !  Todo: move n_msh_pts_p3 later
      ! call list_node_P3 (n_msh_el, n_msh_pts, nodes_per_el, n_edge, n_msh_pts_p3, table_nod, &
      !    table_N_E_F,  visited)

      ! !  TODO: what is signif of this quanitty?
      ! n_ddl = n_edge + n_face + n_msh_pts_p3


      ! if (debug .eq. 1) then
      !    write(ui_out,*) "py_calc_modes.f: n_msh_pts, n_msh_el = ", n_msh_pts, n_msh_el
      !    write(ui_out,*) "py_calc_modes.f: n_msh_pts_p3 = ", n_msh_pts_p3
      !    write(ui_out,*) "py_calc_modes.f: n_vertex, n_edge, n_face,", " n_msh_el = ", &
      !       (n_msh_pts - n_edge), n_edge, n_face, n_msh_el
      !    write(ui_out,*) "py_calc_modes.f: 2D case of the Euler &
      !    & characteristic: V-E+F=1-(number of holes)"
      !    write(ui_out,*) "py_calc_modes.f: Euler characteristic: V - E + F &
      !    &= ", (n_msh_pts - n_edge) - n_edge + n_face
      ! endif


      ! !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
      ! !  Should be using c_dwork for x_E_F ?
      ! call type_node_edge_face (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, type_nod, table_nod, &
      !    table_N_E_F, visited , type_N_E_F, mesh_xy, xy_N_E_F )


      ! !  Fills: type_N_E_F(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
      ! call get_coord_p3 (n_msh_el, n_msh_pts, nodes_per_el, n_ddl, table_nod, type_nod, &
      !    table_N_E_F, type_N_E_F, mesh_xy, xy_N_E_F, visited)



      ! ! From this point ip_visited is unused.

      ! deallocate(visited)

      call set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el, mesh_xy, nodes_per_el, &
         type_nod, table_nod, n_ddl, neq,  xy_N_E_F,  &
         type_N_E_F, m_eqs, debug, &
         iperiod_N, iperiod_N_E_F, inperiod_N, inperiod_N_E_F)

      ! We no longer need type_N_E_F, could deallocate

      !Now we know neq


      !  Sparse matrix CSR setup



         call make_csr_arrays(n_msh_el, n_ddl, neq, table_N_E_F, &
            m_eqs, nonz, v_row_ind, v_col_ptr, debug, errco, emsg);
         RETONERROR(errco)


      !  ----------------------------------------------------------------
      !  convert from 1-based to 0-based
      !  ----------------------------------------------------------------
      !  The CSC indexing, i.e., ip_col_ptr, is 1-based
      !  (but valpr.f will change the CSC indexing to 0-based indexing)


      v_row_ind = v_row_ind - 1
      v_col_ptr = v_col_ptr - 1

      i_base = 0


      write(ui_out,*)
      write(ui_out,*) "-----------------------------------------------"


      vacwavenum_k0 = 2.0d0*D_PI/lambda


      call  check_materials_and_fem_formulation(E_H_field,n_elt_mats, &
         vacwavenum_k0, v_refindex_n, eps_eff, n_core, pp, qq, debug, ui_out, errco, emsg)
      RETONERROR(errco)


      !  Main eigensolver
      write(ui_out,*) "EM FEM: "

      !  Assemble the coefficient matrix A and the right-hand side F of the
      !  finite element equations

      write(ui_out,'(A,A)') "   - assembling linear system:"
      call clock_spare%reset()


      ! These had to wait till we new nonz
      call complex_alloc_1d(mOp_stiff, nonz, 'mOp_stiff', errco, emsg); RETONERROR(errco)
      call complex_alloc_1d(mOp_mass, nonz, 'mOp_mass', errco, emsg); RETONERROR(errco)

      write(*,*) 'write to lmat', nonz
      mOp_stiff(1)= C_ZERO
      mOp_stiff(nonz)= C_ZERO
      write(*,*) 'write to lmat'



      !  Build the actual matrices A (mOp_stiff) and M(mOp_mass) for the arpack solving.

      call asmbly (bdy_cdn, i_base, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, &
         shift_ksqr, bloch_vec, n_elt_mats, pp, qq, &
         table_nod, table_N_E_F, type_el, &
         m_eqs, iperiod_N, iperiod_N_E_F, &
         mesh_xy,  xy_N_E_F,  nonz,  &
         v_row_ind, v_col_ptr, &
         mOp_stiff, mOp_mass )


      dim_krylov = 2*n_modes + n_modes/2 +3


      write(ui_out,'(A,i9,A)') '      ', n_msh_el, ' mesh elements'
      write(ui_out,'(A,i9,A)') '      ', n_msh_pts, ' mesh nodes'
      write(ui_out,'(A,i9,A)') '      ', neq, ' linear equations'
      write(ui_out,'(A,i9,A)') '      ', nonz, ' nonzero elements'
      write(ui_out,'(A,f9.3,A)') '      ', nonz/(1.d0*neq*neq)*100.d0, ' % sparsity'
      write(ui_out,'(A,i9,A)') '      ', neq*(dim_krylov+6)*16/2**20, ' MB est. working memory '

      write(ui_out,'(/,A,A)') '       ', clock_spare%to_string()





      !  This is the main solver.
      !  On completion:
      !  unshifted unsorted eigenvalues are in v_evals_beta[1..n_modes]
      !  eigvectors are in arp arp_evecs



      write(ui_out,'(/,A)') "  - solving linear system: "




      call complex_alloc_2d(arp_evecs, neq, n_modes, 'arp_evecs', errco, emsg); RETONERROR(errco)


      write(ui_out,'(/,A)') "      solving eigensystem"
      call clock_spare%reset()

      call valpr_64( &
         i_base, dim_krylov, n_modes, neq, itermax,  &
         arp_tol, nonz, &
         debug, errco, emsg, &
         v_row_ind, v_col_ptr, &
         mOp_stiff, mOp_mass, &
         v_evals_beta, arp_evecs)
      RETONERROR(errco)





      write(ui_out,'(A,A)') '         ', clock_spare%to_string()


      write(ui_out,'(/,A)') "      assembling modes"
      call clock_spare%reset()

      call rescale_and_sort_eigensolutions(n_modes, shift_ksqr, v_evals_beta, v_eig_index)


      !  The eigenvectors will be stored in the array sol
      !  The eigenvalues and eigenvectors are renumbered
      !  using the permutation vector v_eig_index
      call array_sol ( bdy_cdn, n_modes, n_msh_el, n_msh_pts, n_ddl, neq, nodes_per_el, &
         n_core, bloch_vec, v_eig_index, table_nod, table_N_E_F, type_el, &
         m_eqs, iperiod_N, iperiod_N_E_F, &
         mesh_xy, xy_N_E_F, v_evals_beta, mode_pol, arp_evecs, &
         m_evecs, errco, emsg)
      RETONERROR(errco)


      !  Calculate energy in each medium (typ_el)
      call mode_energy (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_core, &
         table_nod, type_el, n_elt_mats, eps_eff,&
         mesh_xy, m_evecs, v_evals_beta, mode_pol)


      !  Doubtful that this check is of any value: delete?
      !  call check_orthogonality_of_em_sol(n_modes, n_msh_el, n_msh_pts, n_elt_mats, pp, table_nod, &
      !  type_el, mesh_xy, v_evals_beta, m_evecs, &!v_evals_beta_pri, m_evecs_pri,
      !  overlap_L, overlap_file, debug, ui_out, &
      !  pair_warning, vacwavenum_k0, errco, emsg)
      !  RETONERROR(errco)


      !  Should this happen _before_ check_ortho?


      !  The z-component must be multiplied by -ii*beta in order to
      !  get the physical, un-normalised z-component
      !  (see Eq. (25) of the JOSAA 2012 paper)
      !  TODO: is this really supposed to be x i beta , or just x beta  ?
      do i_md=1,n_modes
         m_evecs(3,:,i_md,:) = C_IM_ONE * v_evals_beta(i_md) * m_evecs(3,:,i_md,:)
      enddo


      call array_material_EM (n_msh_el, n_elt_mats, v_refindex_n, type_el, ls_material)

      !  Normalisation. Can't use this if we don't do check_ortho.  Not needed
      !  call normalise_fields(n_modes, n_msh_el, nodes_per_el, m_evecs, overlap_L)


      !if (debug .eq. 1) then
      !  write(ui_out,*) "py_calc_modes.f: CPU time for normalisation :", (time2_J-time1_J)
      !endif
      !
      !  Orthonormal integral
      !  if (debug .eq. 1) then
      !  write(ui_out,*) "py_calc_modes.f: Product of normalised field"
      !  overlap_file = "Orthogonal_n.txt"
      !  call get_clocks( systime1_J, time1_J)
      !  call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_elt_mats, pp, table_nod, &
      !  type_el, mesh_xy, v_evals_beta, v_evals_beta_pri, m_evecs, m_evecs_pri, overlap_L, overlap_file, debug, &
      !  pair_warning, vacwavenum_k0)
      !  call get_clocks( systime2_J, time2_J)
      !  write(ui_out,*) "py_calc_modes.f: CPU time for orthogonal :", (time2_J-time1_J)
      !  endif
      !

      deallocate(v_eig_index, xy_N_E_F, overlap_L, arp_evecs)
      deallocate(mOp_stiff, mOp_mass)
      deallocate(v_row_ind, v_col_ptr)


      write(ui_out,'(A,A)') '         ', clock_spare%to_string()

      write(ui_out,*) "-----------------------------------------------"

      !  call report_results_em(debug, ui_out, &
      !  n_msh_pts, n_msh_el, &
      !  time1, time2, time_fact, time_arpack,  time1_postp, time2_postp, &
      !  lambda, e_h_field, bloch_vec, bdy_cdn,  &
      !  int_max, cmplx_max, cmplx_used,  n_core, n_conv, n_modes, &
      !  n_elt_mats, neq, dim_krylov, &
      !  shift_ksqr, v_evals_beta, eps_eff, v_refindex_n)



   end subroutine calc_em_modes_impl








   subroutine check_materials_and_fem_formulation(E_H_field,n_elt_mats, &
      vacwavenum_k0, v_refindex_n, eps_eff, n_core, pp, qq, debug, ui_out, errco, emsg)

      integer(8), intent(in) :: E_H_field, debug
      integer(8), intent(in) :: n_elt_mats, ui_out
      double precision, intent(in):: vacwavenum_k0
      complex(8), intent(in) :: v_refindex_n(n_elt_mats)
      complex(8), intent(out) :: eps_eff(n_elt_mats)

      integer(8), intent(out) :: n_core(2)
      complex(8), intent(out) :: pp(n_elt_mats), qq(n_elt_mats)
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
      do i=1,n_elt_mats-1

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

   subroutine check_orthogonality_of_em_sol(n_modes, n_msh_el, n_msh_pts, n_elt_mats, pp, table_nod, &
      type_el, mesh_xy, v_evals_beta, m_evecs, &
   !v_evals_beta_pri, m_evecs_pri, &
      overlap_L, overlap_file, debug, ui_out, pair_warning, vacwavenum_k0, errco, emsg)

      use numbatmod
      logical pair_warning

      integer(8), intent(in) :: n_modes, debug, ui_out
      integer(8), intent(in) :: n_msh_pts,  n_msh_el, n_elt_mats
      complex(8) pp(n_elt_mats)

      integer(8), intent(out) :: table_nod(nodes_per_el, n_msh_el)
      integer(8), intent(out) :: type_el(n_msh_el)
      double precision, intent(out) :: mesh_xy(2,n_msh_pts)
      double precision vacwavenum_k0

      complex(8), target, intent(out) :: v_evals_beta(n_modes)
      complex(8), target, intent(out) :: m_evecs(3,nodes_per_el+7,n_modes,n_msh_el)

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

      call orthogonal (n_modes, n_msh_el, n_msh_pts, nodes_per_el, n_elt_mats, pp, table_nod, &
         type_el, mesh_xy, v_evals_beta, m_evecs, &
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
      n_elt_mats, neq, dim_krylov, &
      shift_ksqr, v_evals_beta, eps_eff, v_refindex_n)



      use numbatmod

      integer(8) debug, ui_out, e_h_field, bdy_cdn
      integer(8) int_max, cmplx_max, cmplx_used, int_used, real_max,  n_msh_pts, n_msh_el
      double precision bloch_vec(2), lambda
      double precision time1, time2, start_time, end_time, time_fact, time_arpack, time1_postp
      integer(8) n_conv, n_modes, n_elt_mats, nonz,  n_core(2), neq, dim_krylov
      character(len=FNAME_LENGTH)  log_file
      complex(8), intent(in) :: shift_ksqr
      complex(8), target, intent(out) :: v_evals_beta(n_modes)
      complex(8) eps_eff(n_elt_mats)

      complex(8), intent(in) ::  v_refindex_n(n_elt_mats)

      complex(8) z_tmp
      integer(8) i

      !  TODO: hook these up if needed
      cmplx_max = 0
      int_max = 0
      int_used =0
      nonz = 0
      cmplx_used = 0
      real_max = 0


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

         write(26,*)
         write(26,*) "n_modes, dim_krylov, n_conv = ", n_modes, dim_krylov, n_conv
         !write(26,*) "nonz, n_msh_pts*n_modes, ", "nonz/(n_msh_pts*n_modes) = ", nonz, &
         !  n_msh_pts*n_modes, dble(nonz)/dble(n_msh_pts*n_modes)

         !  write(26,*) "len_skyl, n_msh_pts*n_modes, len_skyl/(n_msh_pts*n_modes) = ",
         !  *   len_skyl, n_msh_pts*n_modes, dble(len_skyl)/dble(n_msh_pts*n_modes)

         write(26,*)
         do i=1,n_modes
            write(26,"(i4,2(g22.14),g18.10)") i, v_evals_beta(i)
         enddo
         write(26,*)
         write(26,*) "n_core = ", n_core
         write(26,*) "eps_eff = ", (eps_eff(i),i=1,n_elt_mats)
         write(26,*) "v_refindex_n = ", (v_refindex_n(i),i=1,n_elt_mats)
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




   subroutine rescale_and_sort_eigensolutions(n_modes, shift_ksqr, v_evals_beta, v_eig_index)

      integer(8), intent(in) :: n_modes
      complex(8), intent(in) :: shift_ksqr
      complex(8) :: v_evals_beta(:)
      integer(8), dimension(:), allocatable :: v_eig_index

      integer(8) i

      complex(8) z_beta

      !TODO: make a function. Turn beta^2 raw eig into actual beta
      do i=1,n_modes
         !  z_tmp0 = v_evals_beta(i)
         !  z_tmp = 1.0d0/z_tmp0+shift_ksqr
         !  z_beta = sqrt(z_tmp)

         z_beta = sqrt(1.0d0/v_evals_beta(i)+shift_ksqr )

         !  Mode classification - we want the forward propagating mode
         if (abs(imag(z_beta)/z_beta) .lt. 1.0d-8) then
            !  re(z_beta) > 0 for forward propagating mode
            if (dble(z_beta) .lt. 0) z_beta = -z_beta
         else
            !  im(z_beta) > 0 for forward decaying evanescent mode  !rarely relevant for us
            if (imag(z_beta) .lt. 0) z_beta = -z_beta
         endif

         v_evals_beta(i) = z_beta
      enddo



      !  order v_evals_beta by magnitudes and store in v_eig_index
      call z_indexx (n_modes, v_evals_beta, v_eig_index)
   end subroutine





end module calc_em_impl




