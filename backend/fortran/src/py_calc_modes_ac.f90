#include "numbat_decl.h"

 !  q_ac :   acoustic wave number (q_ac)
 !  n_modes:  desired number of solved acoustic modes
 !  n_msh_pts:  number of nodes in mesh
 !  n_msh_el:   number of (triang) elements in mesh
 !  n_v_el_material:  number of types of material
 !  v_nd_physindex:    ??
 !  elnd_to_mshpt:
 !  v_el_material:
 !  v_nd_xy
 !  v_eigs_nu:  eigen frequencies nu=omega/(2D_PI) for each mode
 !  femsol_ac:
 !  poln_fracs:


module calc_ac_impl

   use numbatmod
   use alloc

   use class_stopwatch
   use class_MeshRaw
   use class_SparseCSC_AC

contains

   subroutine calc_ac_modes_impl(n_modes, q_ac, dimscale_in_m, shift_nu, &
      bdy_cdn, itermax, tol, debug,  &
      symmetry_flag,  c_tensor, rho, build_acmesh_from_emmesh, &
      mesh_file, n_msh_pts, n_msh_el,n_elt_mats, &
      elnd_to_mshpt, v_el_material, v_nd_physindex,  v_nd_xy, &
      v_eigs_nu, femsol_ac, poln_fracs, nberr)



      integer(8), intent(in) :: n_modes

      complex(8), intent(in) :: q_ac
      double precision, intent(in) :: dimscale_in_m
      complex(8), intent(in) :: shift_nu
      integer(8), intent(in) :: bdy_cdn, itermax, debug
      double precision, intent(in) :: tol
      integer(8), intent(in) :: symmetry_flag, build_acmesh_from_emmesh
      integer(8), intent(in) :: n_elt_mats

      complex(8), intent(in) :: c_tensor(6,6,n_elt_mats)
      complex(8), intent(in) :: rho(n_elt_mats)

      character(len=FNAME_LENGTH), intent(in)  :: mesh_file
      integer(8), intent(in) :: n_msh_pts, n_msh_el

      integer(8), intent(in) :: v_nd_physindex(n_msh_pts)

      integer(8), intent(inout) :: v_el_material(n_msh_el)
      integer(8), intent(inout) :: elnd_to_mshpt(P2_NODES_PER_EL, n_msh_el)

      double precision, intent(inout) ::  v_nd_xy(2,n_msh_pts)

      complex(8), intent(out), target :: v_eigs_nu(n_modes)
      complex(8), intent(out), target :: femsol_ac(3,P2_NODES_PER_EL,n_modes,n_msh_el)
      complex(8), intent(out) :: poln_fracs(4,n_modes)


      type(NBError) nberr

      ! locals

      type(MeshRawAC) mesh_raw
      type(MeshEntitiesAC) entities
      type(SparseCSC_AC) cscmat


      integer(8) :: errco
      character(len=EMSG_LENGTH) :: emsg


      complex(8), dimension(:,:), allocatable :: arpack_evecs


      integer(8) i
      integer(8) ui_out,  namelength


      double precision dim_x, dim_y

      complex(8) shift_omsq
      integer(8)  i_base


      !  Variable used by valpr
      integer(8) ltrav, n_conv
      complex(8) z_beta, z_tmp, z_tmp0
      integer(8), dimension(:), allocatable :: v_eig_index


      !  Variable used by valpr
      integer(8)  nvect

      !  Names and Controls

      character(len=FNAME_LENGTH)  gmsh_file, log_file, gmsh_file_pos


      type(Stopwatch) :: clock_main, clock_spare

      integer(8) :: is_em



      errco = 0
      emsg = ""

      ui_out = stdout


      !  nvect = 2*n_modes + n_modes/2 +3
      nvect = 3*n_modes + 3


      call mesh_raw%allocate(n_msh_pts, n_msh_el, n_elt_mats, nberr)
      RET_ON_NBERR(nberr)

      call entities%allocate(n_msh_el, nberr)
      RET_ON_NBERR(nberr)


      call integer_nalloc_1d(v_eig_index, n_modes, 'v_eig_index', nberr); RET_ON_NBERR(nberr)

      is_em = 0

      !  clean mesh_format
      namelength = len_trim(mesh_file)
      gmsh_file = mesh_file(1:namelength-5)//'.msh'
      gmsh_file_pos = mesh_file(1:namelength)
      log_file = mesh_file(1:namelength-5)//'-AC.log'
      if (debug .eq. 1) then
         write(*,*) "mesh_file = ", mesh_file
         write(*,*) "gmsh_file = ", gmsh_file
      endif


      dim_x = dimscale_in_m
      dim_y = dimscale_in_m
      shift_omsq= (2*D_PI*shift_nu)**2



      call clock_main%reset()


      if (build_acmesh_from_emmesh .eq. 0) then  ! NEVER HAPPENS
         call construct_fem_node_tables_ac (mesh_file, dim_x, dim_y, n_msh_el, n_msh_pts, &
            P2_NODES_PER_EL, n_elt_mats, v_nd_xy, v_nd_physindex, v_el_material, elnd_to_mshpt, errco, emsg)
         call nberr%set(errco, emsg); RET_ON_NBERR(nberr)

      endif


      !  Fills:  MeshRawEM: v_nd_xy, v_nd_physindex, v_el_material, elnd_to_mshpt
      ! This knows the position and material of each elt and mesh point but not their connectedness or edge/face nature

      if (build_acmesh_from_emmesh .eq. 0) then  ! NEVER HAPPENS

         call mesh_raw%construct_node_tables_from_scratch(mesh_file, dimscale_in_m, nberr);
         RET_ON_NBERR(nberr)
      else
         call mesh_raw%construct_node_tables_from_py(v_nd_xy, v_nd_physindex, &
            v_el_material, elnd_to_mshpt, nberr);
         RET_ON_NBERR(nberr)
      endif


      call cscmat%set_boundary_conditions(bdy_cdn, mesh_raw, nberr)
      RET_ON_NBERR(nberr)


      call cscmat%make_csc_arrays(mesh_raw, entities, nberr)
      RET_ON_NBERR(nberr)



      ltrav = 3*nvect*(nvect+2)


      !  The CSC v_eig_indexing, i.e., ip_col_ptr, is 1-based
      !  (but valpr.f will change the CSC v_eig_indexing to 0-based v_eig_indexing)
      i_base = 0

      !#####################  End FEM PRE-PROCESSING  #########################

      write(ui_out,*)
      write(ui_out,*) "-----------------------------------------------"

      write(ui_out,*) "AC FEM: "

      !  Assemble the coefficient matrix K and M of the finite element equations

      write(ui_out,'(A,A)') "   - assembling linear system:"
      call clock_spare%reset()

      call assembly_ac (i_base, shift_omsq, q_ac, rho, c_tensor, &
         mesh_raw, cscmat, symmetry_flag, nberr)
      RET_ON_NBERR(nberr)


      write(ui_out,'(A,i9,A)') '      ', n_msh_el, ' mesh elements'
      write(ui_out,'(A,i9,A)') '      ', n_msh_pts, ' mesh nodes'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_dof, ' linear equations'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_nonz, ' nonzero elements'
      write(ui_out,'(A,f9.3,A)') '      ', cscmat%n_nonz/(1.d0*cscmat%n_dof*cscmat%n_dof)*100.d0, ' % sparsity'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_dof*(nvect+6)*16/2**20, ' MB est. working memory '

      write(ui_out,'(/,A,A)') '       ', clock_spare%to_string()


      write(ui_out,'(/,A)') "  - solving linear system: "
      write(ui_out,'(/,A)') "      solving eigensystem"
      call clock_spare%reset()

      call complex_nalloc_2d(arpack_evecs, cscmat%n_dof, n_modes, 'arpack_evecs', nberr); RET_ON_NBERR(nberr)

      call valpr_64_AC (i_base, nvect, n_modes, cscmat,  itermax, ltrav, tol, &
         n_conv, v_eigs_nu, arpack_evecs, nberr)

      RET_ON_NBERR(nberr)




      if (n_conv .ne. n_modes) then
         write(emsg, '(A,I5,I5)') &
            "py_calc_modes_AC: convergence problem " // &
            "in valpr_64: n_conv != n_modes  ", n_conv, n_modes
         errco = -7
         call nberr%set(errco, emsg); RET_ON_NBERR(nberr)

      endif

      write(ui_out,'(A,A)') '         ', clock_spare%to_string()


      write(ui_out,'(/,A)') "      assembling modes"
      call clock_spare%reset()


      do i=1,n_modes
         z_tmp0 = v_eigs_nu(i)
         z_tmp = 1.0d0/z_tmp0+shift_omsq
         z_beta = sqrt(z_tmp) / (2.0d0 * D_PI)
         !  Frequency (z_beta) should always be positive.
         if (dble(z_beta) .lt. 0) z_beta = -z_beta
         v_eigs_nu(i) = z_beta
      enddo

      call z_indexx_AC (n_modes, v_eigs_nu, v_eig_index)
      !
      !  The eigenvectors will be stored in the array femsol_ac
      !  The eigevalues and eigenvectors will be renumbered
      !  using the permutation vector v_eig_index

      call array_sol_AC (mesh_raw, cscmat, n_modes,   &
         v_eig_index, v_eigs_nu, arpack_evecs, poln_fracs, femsol_ac)




      write(ui_out,'(A,A)') '         ', clock_spare%to_string()
      write(ui_out,*) "-----------------------------------------------"
      write(ui_out,*)

      deallocate(v_eig_index)

   end subroutine calc_ac_modes_impl

end module calc_ac_impl
