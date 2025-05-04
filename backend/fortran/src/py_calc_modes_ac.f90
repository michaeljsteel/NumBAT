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
 !  v_evals_nu:  eigen frequencies nu=omega/(2D_PI) for each mode
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
      bdy_cdn, itermax, arp_tol, &
      symmetry_flag,  c_tensor, rho, build_acmesh_from_emmesh, &
      mesh_file, n_msh_pts, n_msh_el,n_elt_mats, &
      elnd_to_mshpt, v_el_material, v_nd_physindex,  v_nd_xy, &
      v_evals_nu, femsol_ac, poln_fracs, nberr)



      integer(8), intent(in) :: n_modes

      complex(8), intent(in) :: q_ac
      double precision, intent(in) :: dimscale_in_m
      complex(8), intent(in) :: shift_nu
      integer(8), intent(in) :: bdy_cdn, itermax
      double precision, intent(in) :: arp_tol
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

      complex(8), intent(out), target :: v_evals_nu(n_modes)
      complex(8), intent(out), target :: femsol_ac(3,P2_NODES_PER_EL,n_modes,n_msh_el)
      complex(8), intent(out) :: poln_fracs(4,n_modes)


      type(NBError) nberr

      ! locals

      type(MeshRawAC) mesh_raw
      type(MeshEntitiesAC) entities
      type(SparseCSC_AC) cscmat


      character(len=EMSG_LENGTH) :: emsg


      complex(8), dimension(:,:), allocatable :: arpack_evecs

      integer(8) ui_out

      complex(8) shift_omsq
      integer(8) csc_index_offset

      integer(8) shortrun

      integer(8)  dim_krylov

      type(Stopwatch) :: clock_main, clock_spare


      emsg = ""

      ui_out = stdout
      call clock_main%reset()

      call mesh_raw%allocate(n_msh_pts, n_msh_el, n_elt_mats, nberr)
      RET_ON_NBERR(nberr)

      call entities%allocate(n_msh_el, nberr)
      RET_ON_NBERR(nberr)

      if (build_acmesh_from_emmesh .eq. 0) then  ! NEVER HAPPENS

         call construct_fem_node_tables_ac (mesh_file, dimscale_in_m,  n_msh_el, n_msh_pts, &
         n_elt_mats, v_nd_xy, v_nd_physindex, v_el_material, elnd_to_mshpt, nberr)
         RET_ON_NBERR(nberr)

         call mesh_raw%construct_node_tables_from_scratch(mesh_file, dimscale_in_m, nberr);
         RET_ON_NBERR(nberr)
      else
         call mesh_raw%construct_node_tables_from_py(v_nd_xy, v_nd_physindex, &
         v_el_material, elnd_to_mshpt);
      endif


      call cscmat%set_boundary_conditions(bdy_cdn, mesh_raw, nberr)
      RET_ON_NBERR(nberr)


      call cscmat%make_csc_arrays(mesh_raw, nberr)
      RET_ON_NBERR(nberr)


      write(ui_out,*)
      write(ui_out,*) "-----------------------------------------------"

      write(ui_out,*) "AC FEM: "

      !  Assemble the coefficient matrix K and M of the finite element equations

      write(ui_out,'(A,A)') "   - assembling linear system:"
      call clock_spare%reset()


      !  Build the actual matrices A (cscmat%mOp_stiff) and M(cscmat%mOp_mass) for the arpack solving.

      !  The CSC v_eig_indexing, i.e., ip_col_ptr, is 1-based
      !  (but valpr.f will change the CSC v_eig_indexing to 0-based v_eig_indexing)
      csc_index_offset = 0

      shift_omsq= (2*D_PI*shift_nu)**2
      call assembly_ac (csc_index_offset, shift_omsq, q_ac, rho, c_tensor, &
      mesh_raw, cscmat, symmetry_flag, nberr)
      RET_ON_NBERR(nberr)

      !  dim_krylov = 2*n_modes + n_modes/2 +3
      dim_krylov = 3*n_modes + 3

      write(ui_out,'(A,i9,A)') '      ', n_msh_el, ' mesh elements'
      write(ui_out,'(A,i9,A)') '      ', n_msh_pts, ' mesh nodes'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_dof, ' linear equations'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_nonz, ' nonzero elements'
      write(ui_out,'(A,f9.3,A)') '      ', cscmat%n_nonz/(1.d0*cscmat%n_dof*cscmat%n_dof)*100.d0, ' % sparsity'
      write(ui_out,'(A,i9,A)') '      ', cscmat%n_dof*(dim_krylov+6)*16/2**20, ' MB est. working memory '

      write(ui_out,'(/,A,A)') '       ', clock_spare%to_string()


      write(ui_out,'(/,A)') "  - solving linear system: "
      write(ui_out,'(/,A)') "      solving eigensystem"
      call clock_spare%reset()

      call complex_nalloc_2d(arpack_evecs, cscmat%n_dof, n_modes, 'arpack_evecs', nberr); RET_ON_NBERR(nberr)


      shortrun=0
      call solve_arpack_problem (csc_index_offset, dim_krylov, n_modes, itermax,  arp_tol, cscmat, v_evals_nu, arpack_evecs, nberr, shortrun)
      RET_ON_NBERR(nberr)

      write(ui_out,'(A,A)') '         ', clock_spare%to_string()
      write(ui_out,'(/,A)') "      assembling modes"
      call clock_spare%reset()



      !  The eigedim_krylovors will be stored in the array femsol_ac
      !  The eigevalues and eigedim_krylovors will be renumbered
      !  by increasing q_ac

      call construct_solution_fields_ac (shift_omsq,  n_modes, &
      mesh_raw, cscmat,   &
         v_evals_nu, arpack_evecs, femsol_ac, poln_fracs,  nberr)
         RET_ON_NBERR(nberr)


      write(ui_out,'(A,A)') '         ', clock_spare%to_string()
      write(ui_out,*) "-----------------------------------------------"
      write(ui_out,*)

   end subroutine calc_ac_modes_impl

end module calc_ac_impl
