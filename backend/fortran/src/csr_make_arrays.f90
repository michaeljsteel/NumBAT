#include "numbat_decl.h"

subroutine make_csc_arrays(mesh_raw, entities, cscmat, debug, errco, emsg)

   use numbatmod
   use alloc
   use nbinterfacesb

   use class_MeshRaw
   use class_SparseCSC

   type(MeshRaw) :: mesh_raw
   type(MeshEntities) :: entities
   type(SparseCSC) :: cscmat


   integer(8) debug

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

   ! ------------------------------------------

   integer(8) nonz_max, max_row_len
   integer(8), dimension(:), allocatable  :: iwork

   ! ------------------------------------------

   cscmat%n_nonz=0

   call integer_alloc_1d(cscmat%v_col_ptr, cscmat%n_dof+1, 'v_col_ptr', errco, emsg); RETONERROR(errco)

   write(*,*) 'that enty loc', entities%v_xy(1,181), entities%v_xy(2,181)
   call csr_make_col_ptr_loose (mesh_raw%n_msh_el, entities%n_entities, cscmat%n_dof, entities%v_tags, &
      cscmat%m_eqs, cscmat%v_col_ptr, nonz_max)


   ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
   ! length of v_row_ind is determined inside csr_length and so allocated there
   call csr_length (mesh_raw%n_msh_el, entities%n_entities, cscmat%n_dof,  entities%v_tags, cscmat%m_eqs, &
      cscmat%v_row_ind, cscmat%v_col_ptr, &
      nonz_max, cscmat%n_nonz, max_row_len, debug, errco, emsg)
   RETONERROR(errco)
   ! v_row_ind is now size

   call integer_alloc_1d(iwork, 3*entities%n_entities, 'iwork', errco, emsg);
   RETONERROR(errco)

   call sort_csr (cscmat%n_dof, cscmat%n_nonz, max_row_len, cscmat%v_row_ind, cscmat%v_col_ptr,  iwork)


end subroutine

