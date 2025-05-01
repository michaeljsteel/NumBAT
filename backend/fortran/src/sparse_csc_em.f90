#include "numbat_decl.h"


! Represents a matrix in compressed column storage CSC format.
! size is determined by:
!   n_ddl \lapprox NDDL_EM_0 x n_elts     (some entities appear in several elts)
!   n_eq  \lapprox 3 x n_ddl              (some entities have less than 3 dof)

! The dense matrix M has dimension N x N
! For the EM problem, N=n_dof
! m_eqs [3,n_ddl] indicates which of the possible 3 x n_ddl dof are active
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
! Tehre is a trick where part of the CSC contructcion is made with a CSR function
! with the row/col vectors passed in reverse.
! The full matrix appeares to be symmetric, so this is ok, and gives a correct CSC format at the conclusion.




module class_SparseCSC_EM

   use numbatmod
   use alloc
   use class_MeshRawEM
   use class_PeriodicBCs

   private

   type, public :: SparseCSC_EM

   integer(8) n_dof
   integer(8) n_nonz

   integer(8), dimension(:,:), allocatable :: m_eqs

   integer(8), dimension(:), allocatable :: v_row_ind
   integer(8), dimension(:), allocatable :: v_col_ptr

   complex(8), dimension(:), allocatable :: mOp_stiff
   complex(8), dimension(:), allocatable :: mOp_mass

   contains

   procedure :: set_boundary_conditions => SparseCSC_EM_set_boundary_conditions
   procedure :: bound_cond_em => SparseCSC_EM_bound_cond_em

   procedure :: make_csc_arrays => SparseCSC_EM_make_csc_arrays

   procedure :: make_col_ptr_provisional => SparseCSC_EM_make_col_ptr_provisional
   procedure :: make_arrays_final => SparseSC_make_arrays_final


   procedure :: cscmat_contains_elt_row_col => SparseCSC_EM_cscmat_contains_elt_row_col
   procedure :: dump_csc_arrays => SparseCSC_EM_dump_csc_arrays

   end type SparseCSC_EM

   contains

#include "sparse_csc_em_impl.f90"

end module class_SparseCSC_EM