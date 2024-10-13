#include "numbat_decl.h"

module class_SparseCSR

   use numbatmod
   use alloc

   private

   type, public :: SparseCSR

   integer(8) neq
   integer(8) nonz

   integer(8), dimension(:,:), allocatable :: m_eqs

   contains

   end type SparseCSR

   contains

#include "sparse_csr_impl.f90"

end module