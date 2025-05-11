

!#include "numbat_decl.h"


subroutine SparseCSC_EM_set_boundary_conditions(this, bdy_cdn, mesh_raw,  entities, pbcs, nberr)


   class(SparseCSC_EM) :: this
   type(PeriodicBCs) :: pbcs

   integer(8) :: bdy_cdn !, n_dof
   integer(8) :: debug

   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities

   type(NBError) nberr


   !locals
   double precision, dimension(2,2) :: lat_vecs

   if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call this%bound_cond_em (bdy_cdn, entities, nberr)

   elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)
      debug=0
      call periodic_lattice_vec (mesh_raw%n_msh_pts, mesh_raw%v_mshpt_xy, lat_vecs, debug)

      call periodic_node(mesh_raw%n_msh_elts, mesh_raw%n_msh_pts, &
         P2_NODES_PER_EL, mesh_raw%v_mshpt_physindex, mesh_raw%v_mshpt_xy, pbcs%iperiod_N, &
         pbcs%inperiod_N, mesh_raw%m_elnd_to_mshpt, lat_vecs)

      call periodic_N_E_F (entities%n_entities, entities%v_ety_props, entities%v_xy, pbcs%iperiod_N_E_F, &
         pbcs%inperiod_N_E_F, lat_vecs)

      call periodic_cond ( bdy_cdn, entities%n_entities, this%n_dof, entities%v_ety_props, &
         pbcs%iperiod_N_E_F, this%m_eqs, debug)

   endif

end subroutine



!  dof_props = 0  => interiour ddl (ddl = Degree Of Freedom)
!  dof_props != 0 => boundary ddl

!  bdy_cdn = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!  bdy_cdn = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!  bdy_cdn = 2 => Periodic boundary condition



!  This subroutine set the boundary condition parameters


!  The provided Fortran subroutine bound_cond sets boundary condition parameters for a finite element mesh, depending on the specified boundary condition type. Here's a detailed step-by-step explanation of the code:

!  Purpose:
!  The subroutine assigns the degrees of freedom (DOFs) for nodes, edges, and elements in a finite element mesh based on specified boundary conditions (Dirichlet or Neumann).

!  Parameters:
!  Input Parameters:

!    bdy_cdn: Specifies the type of boundary condition.
!      0 for Dirichlet boundary condition (electric wall condition).
!      1 for Neumann boundary condition (magnetic wall condition).
!      2 for Periodic boundary condition.
!
!  n_ddl: Total number of entities.

!
!  dof_props: An array where:
!  The first row indicates if a DOF is an interior or boundary DOF.
!  The second row indicates the dimensional type (0 for nodes, 1 for edges, 2 for elements).

!  Output Parameters:

!    n_dof: Total number of equations (DOFs).
!    m_eqs: An array mapping each DOF to its equation number, considering the boundary conditions.

!  Local Variables:
!  i: Loop index.
!  i_boundary: Indicates if the current DOF is on the boundary.
!  i_dim: Dimensional type of the current DOF (node, edge, or element).
!
!  Description:
!    Dirichlet Boundary Condition (bdy_cdn = 0):
!      All points (nodes, edges, elements) have degrees of freedom unless they are boundary points.
!      For elements (i_dim = 2): Each element has 3 interior DOFs.
!      For edges (i_dim = 1): Each edge has 3 DOFs unless it is on the boundary.
!      For nodes (i_dim = 0): Each node has 1 DOF unless it is on the boundary.
!
!  Neumann Boundary Condition (bdy_cdn = 1):
!     All points (nodes, edges, elements) have degrees of freedom.
!     For elements (i_dim = 2) and edges (i_dim = 1): Each has 3 DOFs.
!     For nodes (i_dim = 0): Each node has 1 DOF.

!  matrix m_eqs is built up based on properties of the entities established
!    in MeshEntities
!
! dof_props holds the is_boundary and dimensionality of each entity

!  m_eqs assigns an index to each degree of freedom, if it exists, for each entity

subroutine SparseCSC_EM_bound_cond_em (this, bdy_cdn, entities, nberr)

   class(SparseCSC_EM) :: this
   type(MeshEntities) :: entities

   integer(8) bdy_cdn, n_dof
   type(NBError) nberr

   character(len=EMSG_LENGTH) :: emsg

   integer(8) i, is_boundary, i_dim

   call integer_alloc_2d(this%m_eqs, 3_8, entities%n_entities, 'm_eqs', nberr); RET_ON_NBERR(nberr)

   n_dof = 0

   if (bdy_cdn .eq. BCS_DIRICHLET) then  !  all points have a degree of freedom

      do i=1,entities%n_entities
         is_boundary = entities%v_ety_props(1,i)
         i_dim = entities%v_ety_props(2,i)

         if (i_dim .eq. 2) then !  each element is associated with 3 interior Degrees Of Freedom (DOF)
            this%m_eqs(1,i) = n_dof + 1
            this%m_eqs(2,i) = n_dof + 2
            this%m_eqs(3,i) = n_dof + 3
            n_dof = n_dof + 3

         elseif (i_dim .eq. 1) then  !  each edge is associated with 3 Degrees Of Freedom (DOF)
            if (is_boundary .eq. 0) then
               this%m_eqs(1,i) = n_dof + 1
               this%m_eqs(2,i) = n_dof + 2
               this%m_eqs(3,i) = n_dof + 3
               n_dof = n_dof + 3
            else                     !  fixed by boundary so no dof
               this%m_eqs(1,i) = 0
               this%m_eqs(2,i) = 0
               this%m_eqs(3,i) = 0
            endif

         elseif (i_dim .eq. 0) then   !  each node is associated with 1 Degree Of Freedom (DOF)
            if (is_boundary .eq. 0) then
               this%m_eqs(1,i) = n_dof + 1
               this%m_eqs(2,i) = 0
               this%m_eqs(3,i) = 0
               n_dof = n_dof + 1
            else
               this%m_eqs(1,i) = 0
               this%m_eqs(2,i) = 0
               this%m_eqs(3,i) = 0
            endif
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
               "bdy_cdn = ", bdy_cdn, "i = ", i
            call nberr%set(NBERR_BAD_BOUNDARY_CONDITION, emsg)
            return
         endif
      enddo

   elseif(bdy_cdn .eq. BCS_NEUMANN) then !all points have a degree of freedom

      do i=1,entities%n_entities
         i_dim = entities%v_ety_props(2,i)
         if (i_dim .eq. 2 .or. i_dim .eq. 1) then !  Each element or edge is associated with 3 Degrees Of Freedom (DOF)
            this%m_eqs(1,i) = n_dof + 1
            this%m_eqs(2,i) = n_dof + 2
            this%m_eqs(3,i) = n_dof + 3
            n_dof = n_dof + 3
         elseif (i_dim .eq. 0) then
            this%m_eqs(1,i) = n_dof + 1
            this%m_eqs(2,i) = 0
            this%m_eqs(3,i) = 0
            n_dof = n_dof + 1
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
               "bdy_cdn = ", bdy_cdn, "i = ", i

            call nberr%set(NBERR_BAD_BOUNDARY_CONDITION, emsg)
            return
         endif
      enddo
   endif

   this%n_dof = n_dof

end



subroutine SparseCSC_EM_make_csc_arrays(this, bdy_cdn, mesh_raw, entities, pbcs, nberr)

   class(SparseCSC_EM) :: this
   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities
   type(PeriodicBCs) :: pbcs

   type(NBError) nberr

   integer(8) bdy_cdn
    ! ------------------------------------------

   integer(8) n_nonz_max, max_row_len

   call this%set_boundary_conditions(bdy_cdn, mesh_raw, entities, pbcs, nberr);
   RET_ON_NBERR(nberr)


   this%n_nonz=0

   call integer_alloc_1d(this%v_col_ptr, this%n_dof+1, 'v_col_ptr', nberr); RET_ON_NBERR(nberr)

   call this%make_col_ptr_provisional (mesh_raw, entities, n_nonz_max)


   ! v_col_ptr now has the right length for CSC and n_nonz_max is an upper bound for the number of n_nonzeros.
   ! Now get the row_indexes.
   call this%make_arrays_final (mesh_raw, entities, n_nonz_max, max_row_len, nberr);
   RET_ON_NBERR(nberr)


   ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
   ! length of v_row_ind is determined inside csr_length and so allocated there
   !call csr_length (mesh_raw%n_msh_elts, entities%n_entities, this%n_dof,  entities%v_tags, this%m_eqs, &
   !this%v_row_ind, this%v_col_ptr, n_nonz_max, this%n_nonz, max_row_len, errco, emsg)
   !RETONERROR(errco)


   ! ! now have valid construction of minimal length csc matrix
   ! do row =1, this%n_dof
   !    do col=1,  this%n_dof
   !       call this%cscmat_contains_elt_row_col(row, col, rc_exists, val)
   !       if (rc_exists .eq. 0) cycle

   !       ! found one, is it symmetric?
   !       call this%cscmat_contains_elt_row_col(col, row, cr_exists, val)

   !       if (rc_exists .ne. 0) then
   !          !write(*,*) 'csc found sym pair at', row, col
   !       else
   !          write(*,*) 'csc found asym elt at', row, col
   !       endif
   !    enddo
   ! enddo

   !call this%dump_csc_arrays()




   ! At this point, the row_indices for a given column are in random order
   ! Now we sort them column by column so the rows appear in ascending order within each column
   ! This is another reverse passing to a CSR routine
   call sort_csc (this%n_dof, this%n_nonz, max_row_len, this%v_row_ind, this%v_col_ptr)

   !call this%dump_csc_arrays()

   ! Now we know how big the data arrays are
   call complex_alloc_1d(this%mOp_stiff, this%n_nonz, 'cscmat%mOp_stiff', nberr); RET_ON_NBERR(nberr)

   call complex_alloc_1d(this%mOp_mass, this%n_nonz, 'cscmat%mOp_mass', nberr); RET_ON_NBERR(nberr)

end subroutine



! subroutine SparseCSC_EM_cscmat_contains_elt_row_col(this, row, col, found, val)

!    class(SparseCSC_EM) :: this
!    integer(8) row, col
!    integer(8) found

!    integer(8) val
!    integer(8) vct_1, vct_2, ri

!    found = 0
!    if (col .lt. 1 .or. col .gt. this%n_dof) then
!       return
!    endif

!    vct_1 = this%v_col_ptr(col)
!    vct_2 = this%v_col_ptr(col+1)
!    do ri = vct_1, vct_2
!       if (this%v_row_ind(ri) .eq. row) then
!          found = 1
!          val = 1
!       endif
!    enddo

! end subroutine


subroutine SparseCSC_EM_dump_csc_arrays(this)
   class(SparseCSC_EM) :: this
   integer(8) i


   write(*,*) 'CSC col ptrs'
   do i=1,this%n_dof+1
      write(*,*) i, this%v_col_ptr(i)
   enddo


   write(*,*) 'CSC row indices'
   do i=1,this%n_nonz
      write(*,*) i, this%v_row_ind(i)
   enddo


end subroutine



subroutine SparseCSC_EM_make_col_ptr_provisional (this, mesh_raw, entities, n_nonz_max)


   class(SparseCSC_EM) :: this
   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities


   integer(8) n_nonz_max
   !integer(8) ety_tags (14,nel)
   !integer(8) m_eqs(3,n_entty), v_col_ptr(n_dof+1)


   integer(8)  k_el, active_dof, tag, i_nd, i_col, locdof
   integer(8) last_ct, this_ct


   this%v_col_ptr = 0


   !  Count references to each dof
   !  This is equivalent to counting how many elements a given entity falls on.

   !  Here v_col_ptr is just a convenient temporary memory holder.
   !  The contents is not related to its actual definition as the coloumn pointer

   do  k_el=1, mesh_raw%n_msh_elts          ! for every ety at every elt

      do i_nd=1,N_ENTITY_PER_EL
         tag = entities%v_tags(i_nd,k_el)      ! find its tag

         do locdof = 1,3                   ! for each of its 3 possible dof,
            active_dof = this%m_eqs(locdof, tag)    ! find the index of that dof, n_nonzero means active

            ! count the number of times the dof is encountered
            if (active_dof .ne. 0) this%v_col_ptr(active_dof) = this%v_col_ptr(active_dof)+1
         enddo

      enddo
   enddo


   ! v_col_ptr now contains the number of mentions of each dof, ie the number of equations it is involved in
   ! This is more or less the number of triangles it lies on,
   ! The largest values occur at elt entities 5,6 or probably 7 which are vertices when many triangles meet
   ! This is usually of order 6-10


   !  Compressed Column Storage (CSC): determine the column pointer
   !   Each v_col_ptr(j)  = 1 + number of n_nonzeros left of column j
   !                      = v_col_ptr(j-1) + number of n_nonzeros column (j-1)
   !                                  where we set the virtual zero column as v_col_ptr(0)=1
   last_ct = this%v_col_ptr(1)
   this%v_col_ptr(1) = 1          !  The first column begins with element 1 (every dof has some self energy so this is always true )
   do i_col=2,this%n_dof+1
      this_ct = this%v_col_ptr(i_col)    ! # of neigbours of entity i_col (including self)

      ! Set v_col_ptr(i_col) by counting the possible elements in the previous column
      ! An upper bound for how many dof this one might interact with is:
      !     3 dof for each N_ENTITY_PER_EL entities in its own elt, + 3 dof for each of the (N_ENTITY_PER_EL-1) neighbour entities
      ! Some of these will be redundant by double counting, for example the edge and vertex nodes on two adjacent triangles
      !   which are themselves touching
      ! This is still _much_ smaller than interacting with every dof in the mesh
      this%v_col_ptr(i_col) = this%v_col_ptr(i_col-1) + 3*N_ENTITY_PER_EL + 3*(N_ENTITY_PER_EL-1)*(last_ct-1)
      last_ct = this_ct
   enddo

   n_nonz_max = this%v_col_ptr(this%n_dof+1) - 1

end subroutine


! subroutine SparseSC_make_arrays_final (this, mesh_raw, entities, n_nonz_max, max_row_len, errco, emsg)


!    class(SparseCSC_EM) :: this
!    type(MeshRawEM) :: mesh_raw
!    type(MeshEntities) :: entities

!    integer(8) n_nonz_max, max_row_len
!    integer(8) errco
!    character(len=EMSG_LENGTH) emsg

!    ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
!    ! length of v_row_ind is determined inside csr_length and so allocated there

!    write(*,*) 'maf 1'

!    call csr_length (mesh_raw%n_msh_elts, entities%n_entities, this%n_dof,  entities%v_tags, this%m_eqs, &
!       this%v_row_ind, this%v_col_ptr, n_nonz_max, this%n_nonz, max_row_len, errco, emsg)
!    RETONERROR(errco)



!end subroutine



! This one is written in CSR format
! row/col names seem backward
! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine SparseSC_make_arrays_final (this, mesh_raw, entities, n_nonz_max, max_row_len, nberr)


   class(SparseCSC_EM) :: this
   type(MeshRawEM) :: mesh_raw
   type(MeshEntities) :: entities

   integer(8) n_nonz_max, max_row_len
   type(NBError) nberr


   ! --------------------------------------------

   integer(8), dimension(:), allocatable :: row_ind_tmp

   integer(8) i, j, j_nd, k, k1,  i_nd, i_tag, j_tag, i_locdof, j_locdof
   integer(8) k_el, i_dof, j_dof
   integer(8) col_start, col_end, row_len
   integer(8) col_start2, col_end2, ui, stored
   integer(8) n_nonz


   ui = stdout


   call integer_alloc_1d(row_ind_tmp, n_nonz_max, 'row_ind_tmp', nberr); RET_ON_NBERR(nberr)

   row_ind_tmp = 0

   ! This code was converted from one in CSR format
   !  Determination of the row indices

   !write(*,*) 'Total dof is ', n_entty, n_dof, n_nonz_max
   n_nonz = 0
   do k_el=1,mesh_raw%n_msh_elts                    ! for each element

      do i_nd=1,N_ENTITY_PER_EL               !   and its 14 entities
         i_tag = entities%v_tags(i_nd, k_el)

         do i_locdof=1,3                     !   and their 3 potential dof
            i_dof = this%m_eqs(i_locdof, i_tag)   !   When n_nonzero, this is the row number for this dof

            if (i_dof .eq. 0) cycle     ! an inactive dof for this entity, go around again



            col_start = this%v_col_ptr(i_dof)          ! range of elt indices which are in this row
            col_end = this%v_col_ptr(i_dof+1) - 1

            !write(*,*) 'looking for partners of', k_el, i_nd, i_tag, i_locdof, &
            !  i_dof, 'in cols', col_start, col_end, col_end-col_start+1

            do j_nd=1,N_ENTITY_PER_EL
               j_tag = entities%v_tags(j_nd,k_el)

               do j_locdof=1,3
                  j_dof = this%m_eqs(j_locdof, j_tag)

                  if (j_dof .eq. 0) cycle

                  !  Store the entry (i_dof,j_dof) if it's not already found
                  stored = 0

                  do k=col_start,col_end
                     if (row_ind_tmp(k) .eq. 0) then  !goto 20 ! if we find a zero, we've seen all of the n_nonzeros, and none were a match
                        !  No entry exists for (i_dof,j_dof); create new one
                        n_nonz = n_nonz + 1

                        row_ind_tmp(k) = j_dof
                        stored = 1
                        exit
                     endif

                     if (row_ind_tmp(k) .eq. j_dof) then !goto 30  ! already stored, bail out
                        stored=1
                        exit
                     endif
                  enddo

                  if (stored .eq. 0) then ! shouldn't have got here
                     call nberr%set(NBERROR_118,  "csr_length: There is a problem with row/col indexing!")
                     return
                  endif

               enddo
            enddo
         enddo
      enddo
   enddo



   ! The full coupling matrix is dimension n_dof x n_dof
   ! At most n_nonz_max of these can be n_nonzero based on the pairs lying on adjacent elements
   ! We now have a valid CSR indexing of this subset of n_nonz_max pairs

   ! But more of these can be eliminated. Not quite sure why.


   !  squeeze away the zero entries
   !  added so as to handle more type of domains/meshes



   do i=1,this%n_dof-1
      col_start = this%v_col_ptr(i)
      col_end = this%v_col_ptr(i+1) - 1

      do j=col_start,col_end
         if(row_ind_tmp(j) .ne. 0) cycle   ! this elt is busy, check the next

         ! we've found a zero element in this col all the rest of this col will be zeros, so remove them all

         col_start2 = this%v_col_ptr(i) + j - col_start
         this%v_col_ptr(i+1) = col_start2         ! bring the start of the next col forward
         col_end2 = this%v_col_ptr(i+2) - 1       ! find the end of that next col
         do k=col_end+1,col_end2
            k1 = col_start2 + k - (col_end+1)  ! shuffle the columns in that col forward into the empty space
            row_ind_tmp(k1) = row_ind_tmp(k)
            row_ind_tmp(k) = 0
         enddo

         exit  ! this row is done, go back to the next

      enddo


   enddo



   ! squeeze the last col
   i = this%n_dof
   col_start = this%v_col_ptr(i)
   col_end = this%v_col_ptr(i+1) - 1
   do j=col_start,col_end
      if(row_ind_tmp(j) .eq. 0) then
         col_start2 = this%v_col_ptr(i) + j - col_start
         this%v_col_ptr(i+1) = col_start2
         exit
      endif
   enddo


   ! Find the longest row
   max_row_len = 0
   do i=1,this%n_dof
      col_start = this%v_col_ptr(i)
      col_end = this%v_col_ptr(i+1) - 1
      row_len = col_end - col_start + 1
      if (row_len .gt. max_row_len) max_row_len = row_len
   enddo



   ! Now we know n_nonz
   call integer_alloc_1d(this%v_row_ind, n_nonz, 'this%v_row_ind', nberr); RET_ON_NBERR(nberr)

   this%v_row_ind(1:n_nonz) = row_ind_tmp(1:n_nonz)


!   deallocate(row_ind_tmp)

   this%n_nonz = n_nonz

end



  !  convert from 1-based to 0-based
   !  ----------------------------------------------------------------
   !  Our CSC indexing so far, i.e., ip_col_ptr, has been 1-based
   !  But external calls to functions like valpr.f will need 0-based indexing.
   !  So the indices need shifting by one.
   !  (But not the data arrays: c[0] and fortran[1] refer to the same location, so that just works)


! TODO: maintain current offset flag to know if we need to apply shift
subroutine SparseCSC_EM_adjust_for_zero_offset_indexing(this)

   class(SparseCSC_EM) :: this

   this%v_row_ind = this%v_row_ind - 1
   this%v_col_ptr = this%v_col_ptr - 1

end