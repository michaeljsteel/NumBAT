#include "numbat_decl.h"


! Represents a matrix in compressed column storage CSC format.
! size is determined by:
!   n_ddl \lapprox NDDL_EM_0 x n_elts     (some entities appear in several elts)
!   n_eq  \lapprox 3 x n_ddl              (some entities have less than 3 dof)

! The dense matrix M has dimension N x N
! For the EM problem, N=n_dof
! m_global_dofs [3,n_ddl] indicates which of the possible 3 x n_ddl dof are active
!       n_dof is the number of nonzero entries, ie. number of dof
!       Its entries assign an index to all active dof which is the column number of its equation

! M has structure
!  eq1_for_dof1: [dof_1  dof_2  dof_3  ...]
!  eq2_for_dof1: [dof_1  dof_2  dof_3  ...]
!  eqN1_for_dof1: [dof_1  dof_2  dof_3  ...]
!  eq1_for_dof2: [dof_1  dof_2  dof_3  ...]
!  eq2_for_dof1: [dof_1  dof_2  dof_3  ...]
!  ...
!  eq1_for_dof3: [dof_1  dof_2  dof_3  ...]
!  ...
!

! We can construct row and column ptr arrays even though we don't yet know the actual values of the
! nonzero entries of M, only where they are.
! M is represented by three arrays  v_M, v_row_ind, v_col_ptr
! v_M       : length nonz, the number of nonzero elements
!              contains the nonzero data in column-major order
! v_row_ind : length nonz
!              contains the row of each entry in v_M
! v_col_ptr : length N+1 = 3 x NNDL_0_EM + 1
!           encodes the index in V where a new column starts
!           ie, v_col_ptr(j=1..N) contains 1 plus the number of nonzero elements left of column j
!               v_col_ptr(N+1)= nonz
!
!
! There is a trick where part of the CSC contsruction is made with a CSR function
! with the row/col vectors passed in reverse.
! The full matrix appears to be symmetric, so this is ok, and gives a correct CSC format at the conclusion.




module class_SparseCSC_AC

   use numbatmod
   use alloc
   use class_Mesh
   use class_PeriodicBCs

   private

   type, public :: SparseCSC_AC

   integer(8) n_dof
   integer(8) n_nonz

   ! maps xyz_locdof (local DOF) and mesh_pt to absolute dof
   ! size: [3, n_msh_pts]
   ! ranges: [0..n_dof]
   ! 0 value means no dof is associated with this site and local dof, due to bdy conditions
   ! natural name would be m_global_dof
   integer(8), dimension(:,:), allocatable :: m_global_dofs

   ! row indices for each nonzero element in the sparse representation
   ! size: [n_nonz], ranges: [1..n_dof]
   integer(8), dimension(:), allocatable :: v_row_ind

   ! starting index for each new column in the sparse representation
   ! size: [<=n_dof+1], ranges: [1..n_dof]
   integer(8), dimension(:), allocatable ::  v_col_ptr

   ! CSC representation of the FEM stiffness K and mass M matrices
   ! (not to be confused with elastic stiffness C)
   ! size: [n_nonz]
   complex(8), dimension(:), allocatable :: mOp_stiff
   complex(8), dimension(:), allocatable :: mOp_mass

   contains

   procedure :: set_boundary_conditions => SparseCSC_AC_set_boundary_conditions
   !procedure :: set_boundary_conditions => SparseCSC_AC_set_bound_cond

   procedure :: make_csc_arrays => SparseCSC_AC_make_csc_arrays

   procedure :: make_col_ptr_provisional => SparseCSC_AC_make_col_ptr_provisional
   procedure :: make_arrays_final => SparseCSC_AC_make_arrays_final

   procedure :: adjust_for_zero_offset_indexing => SparseCSC_AC_adjust_for_zero_offset_indexing


   !procedure :: cscmat_contains_elt_row_col => SparseCSC_AC_cscmat_contains_elt_row_col
   procedure :: dump_csc_arrays => SparseCSC_AC_dump_csc_arrays

   end type SparseCSC_AC

   contains

#include "sparse_csc_ac_impl.f90"

end module class_SparseCSC_AC


