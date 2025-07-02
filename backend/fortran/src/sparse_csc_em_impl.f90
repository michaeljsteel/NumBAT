

subroutine SparseCSC_EM_set_boundary_conditions(this, bdy_cdn, mesh, entities, pbcs, nberr)


   class(SparseCSC_EM) :: this

   integer(8) :: bdy_cdn
   integer(8) :: debug

   type(MeshEM) :: mesh
   type(MeshEntities) :: entities
   type(PeriodicBCs) :: pbcs

   type(NBError) nberr


   !locals
   double precision, dimension(2,2) :: lat_vecs

   if (bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call this%bound_cond_em (bdy_cdn, entities, nberr)

   elseif(bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)
      debug=0
      call periodic_lattice_vec (mesh%n_msh_pts, mesh%v_mshpt_xy, lat_vecs, debug)

      call periodic_node(mesh%n_msh_elts, mesh%n_msh_pts, &
         P2_NODES_PER_EL, mesh%v_mshpt_physindex, mesh%v_mshpt_xy, pbcs%iperiod_N, &
         pbcs%inperiod_N, mesh%m_elnd_to_mshpt, lat_vecs)

      call periodic_N_E_F (entities%n_entities, entities%v_ety_props, entities%v_xy, pbcs%iperiod_N_E_F, &
         pbcs%inperiod_N_E_F, lat_vecs)

      call periodic_cond ( bdy_cdn, entities%n_entities, this%n_dof, entities%v_ety_props, &
         pbcs%iperiod_N_E_F, this%m_global_dofs, debug)

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
!    m_global_dofs: An array mapping each DOF to its equation number, considering the boundary conditions.

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

!  matrix m_global_dofs is built up based on properties of the entities established
!    in MeshEntities
!
! dof_props holds the is_boundary and dimensionality of each entity

!  m_global_dofs assigns an index to each degree of freedom, if it exists, for each entity

! Total DOF:
!   face 3
!   edges 3 x3
!   P3 nodes 10
!   Total: 22  !

subroutine SparseCSC_EM_bound_cond_em (this, bdy_cdn, entities, nberr)

   class(SparseCSC_EM) :: this

   type(MeshEntities) :: entities
   integer(8) bdy_cdn
   type(NBError) nberr

   character(len=EMSG_LENGTH) :: emsg

   integer(8) ety_i, is_boundary, i_dim, n_dof

   call integer_alloc_2d(this%m_global_dofs, 3_8, entities%n_entities, 'm_global_dofs', nberr); 

   RET_ON_NBERR(nberr) 
   n_dof = 0

   if (bdy_cdn .eq. BCS_DIRICHLET) then  !  all points have a degree of freedom

      do ety_i=1,entities%n_entities
         is_boundary = entities%v_ety_props(ETY_PROP_PHYSTYPE, ety_i)
         i_dim = entities%v_ety_props(ETY_PROP_DIMENSION, ety_i)    ! this is not really a dimension, but an identifier: face, P2 edge, or P3 node

         if (i_dim .eq. 2) then ! is a face       ! A face has a set of 3 vec states, like a P2 edge
            !  each ety has 3 dof:
            this%m_global_dofs(1, ety_i) = n_dof + 1
            this%m_global_dofs(2, ety_i) = n_dof + 2
            this%m_global_dofs(3, ety_i) = n_dof + 3
            n_dof = n_dof + 3

         elseif (i_dim .eq. 1) then  ! is a P2 edge, for 3 vec dofs each, making 12 P2-related DOF?
            if (is_boundary .eq. 0) then
               this%m_global_dofs(1, ety_i) = n_dof + 1
               this%m_global_dofs(2, ety_i) = n_dof + 2
               this%m_global_dofs(3, ety_i) = n_dof + 3
               n_dof = n_dof + 3
            else                     !  fixed by boundary so no dof
               this%m_global_dofs(1, ety_i) = 0
               this%m_global_dofs(2, ety_i) = 0
               this%m_global_dofs(3, ety_i) = 0
            endif

         elseif (i_dim .eq. 0) then   !  is a P3 node, 1 dof for each of 10 nodes (since scalar field)
            if (is_boundary .eq. 0) then
               this%m_global_dofs(1, ety_i) = n_dof + 1
               this%m_global_dofs(2, ety_i) = 0
               this%m_global_dofs(3, ety_i) = 0
               n_dof = n_dof + 1
            else
               this%m_global_dofs(1, ety_i) = 0
               this%m_global_dofs(2, ety_i) = 0
               this%m_global_dofs(3, ety_i) = 0
            endif
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
               "bdy_cdn = ", bdy_cdn, "i = ", ety_i
            call nberr%set(NBERR_BAD_BOUNDARY_CONDITION, emsg)
            return
         endif
      enddo

   elseif(bdy_cdn .eq. BCS_NEUMANN) then !all points have a degree of freedom

      do ety_i=1,entities%n_entities
         i_dim = entities%v_ety_props(ETY_PROP_DIMENSION, ety_i)
         if (i_dim .eq. 2 .or. i_dim .eq. 1) then !  Each element or edge is associated with 3 Degrees Of Freedom (DOF)
            this%m_global_dofs(1, ety_i) = n_dof + 1
            this%m_global_dofs(2, ety_i) = n_dof + 2
            this%m_global_dofs(3, ety_i) = n_dof + 3
            n_dof = n_dof + 3
         elseif (i_dim .eq. 0) then
            this%m_global_dofs(1, ety_i) = n_dof + 1
            this%m_global_dofs(2, ety_i) = 0
            this%m_global_dofs(3, ety_i) = 0
            n_dof = n_dof + 1
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
               "bdy_cdn = ", bdy_cdn, "i = ", ety_i

            call nberr%set(NBERR_BAD_BOUNDARY_CONDITION, emsg)
            return
         endif
      enddo
   endif

   this%n_dof = n_dof

end



subroutine SparseCSC_EM_make_csc_arrays(this, bdy_cdn, mesh, entities, pbcs, nberr)

   class(SparseCSC_EM) :: this
   type(MeshEM) :: mesh
   type(MeshEntities) :: entities
   type(PeriodicBCs) :: pbcs

   type(NBError) nberr

   integer(8) bdy_cdn
   ! ------------------------------------------

   integer(8) n_nonz_max, max_col_len

   call this%set_boundary_conditions(bdy_cdn, mesh, entities, pbcs, nberr);
   RET_ON_NBERR(nberr)


   this%n_nonz=0

   call integer_alloc_1d(this%v_col_ptr, this%n_dof+1, 'v_col_ptr', nberr); RET_ON_NBERR(nberr)

   call this%make_col_ptr_provisional (mesh, entities, n_nonz_max)


   ! v_col_ptr now has the right length for CSC and n_nonz_max is an upper bound for the number of n_nonzeros.
   ! Now get the row_indexes.
   call this%make_arrays_final (mesh, entities, n_nonz_max, max_col_len, nberr);
   RET_ON_NBERR(nberr)


   ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
   ! length of v_row_ind is determined inside csr_length and so allocated there
   !call csr_length (mesh%n_msh_elts, entities%n_entities, this%n_dof,  entities%v_tags, this%m_global_dofs, &
   !this%v_row_ind, this%v_col_ptr, n_nonz_max, this%n_nonz, max_col_len, errco, emsg)
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
   call sort_csc (this%n_dof, this%n_nonz, max_col_len, this%v_row_ind, this%v_col_ptr)

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



subroutine SparseCSC_EM_make_col_ptr_provisional (this, mesh, entities, n_nonz_max)


   class(SparseCSC_EM) :: this
   type(MeshEM) :: mesh
   type(MeshEntities) :: entities


   integer(8) n_nonz_max
   !integer(8) ety_tags (14,nel)
   !integer(8) m_global_dofs(3,n_entty), v_col_ptr(n_dof+1)


   integer(8)  el_i, dof, tag, nd_i, i_dof, locdof
   integer(8) last_ct, this_ct


   this%v_col_ptr = 0


   ! Each dof can play multiple roles according to how many elements it participates in
   ! (eg corner>=3, edge=1 or 2)
   ! Counting these gives an upper bound to the nonzero elements of the FEM matrices

   do  el_i=1, mesh%n_msh_elts   ! for every element (triangle)

      do nd_i=1,N_ENTITY_PER_EL      ! and all nodes/etys in the elt

         ! find global mesh point tag
         tag = entities%v_tags(nd_i,el_i)

         do locdof = 1,3  ! for each of its 3 possible dof,
            ! find the index of that dof, n_nonzero means active
            dof = this%m_global_dofs(locdof, tag)

            ! increment number of roles this dof plays in multiple elements
            if (dof .ne. 0) this%v_col_ptr(dof) = this%v_col_ptr(dof)+1
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

   ! Set v_col_ptr(i_dof) by counting the possible elements in the previous column
   ! An upper bound for how many dof this one might interact with is:
   !     3 dof for each N_ENTITY_PER_EL entities in its own elt, + 3 dof for each of the (N_ENTITY_PER_EL-1) neighbour entities
   ! Some of these will be redundant by double counting, for example the edge and vertex nodes on two adjacent triangles
   !   which are themselves touching
   ! This is still _much_ smaller than interacting with every dof in the mesh

   last_ct = this%v_col_ptr(1)
   this%v_col_ptr(1) = 1          !  The first column begins with element 1 (every dof has some self energy so this is always true )
   do i_dof=2,this%n_dof+1
      this_ct = this%v_col_ptr(i_dof)    ! # of neigbours of entity i_dof (including self)

      this%v_col_ptr(i_dof) = this%v_col_ptr(i_dof-1) + 3*N_ENTITY_PER_EL + 3*(N_ENTITY_PER_EL-1)*(last_ct-1)
      last_ct = this_ct
   enddo

   n_nonz_max = this%v_col_ptr(this%n_dof+1) - 1

end subroutine



! This one is written in CSR format
! row/col names seem backward
! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine SparseSC_make_arrays_final (this, mesh, entities, n_nonz_max, max_col_len, nberr)


   class(SparseCSC_EM) :: this
   type(MeshEM) :: mesh
   type(MeshEntities) :: entities

   integer(8) n_nonz_max, max_col_len
   type(NBError) nberr


   ! --------------------------------------------

   integer(8), dimension(:), allocatable :: row_ind_tmp
   integer(8) i, j, j_nd, k, k1,  nd_i, i_tag, j_tag, i_locdof, j_locdof
   integer(8) el_i, i_dof, j_dof
   integer(8) row_lo, row_hi, row_len
   integer(8) row_lo2, row_hi2, ui
   integer(8) n_nonz


   ui = stdout


   call integer_alloc_1d(row_ind_tmp, n_nonz_max, 'row_ind_tmp', nberr); RET_ON_NBERR(nberr)

   row_ind_tmp = 0

   ! This code was converted from one in CSR format
   !  Determination of the row indices

   ! For two dof to interact and require nonzero element in the FEM matrices,
   ! they must both exist on the same element

   n_nonz = 0
   do el_i=1,mesh%n_msh_elts        ! for each element

      do nd_i=1,N_ENTITY_PER_EL         !   and its 14 entities
         i_tag = entities%v_tags(nd_i, el_i)

         do i_locdof=1,3        !   and their 3 potential dof
            i_dof = this%m_global_dofs(i_locdof, i_tag)   !   When n_nonzero, this is the row number for this dof

            if (i_dof .eq. 0) cycle   ! an inactive dof



            ! range of elt indices which are in this col
            row_lo = this%v_col_ptr(i_dof)
            row_hi = this%v_col_ptr(i_dof+1) - 1


            do j_nd=1,N_ENTITY_PER_EL
               j_tag = entities%v_tags(j_nd,el_i)

               do j_locdof=1,3
                  j_dof = this%m_global_dofs(j_locdof, j_tag)

                  if (j_dof .eq. 0) cycle


                  ! Now we know (i_dof, j_dof) can interact and should have a slot
                  ! in the column for i_dof

                  !  Store the entry (i_dof,j_dof) if it's not already found
                  call store_dof_in_csc_row_index(row_lo, row_hi, j_dof, row_ind_tmp, n_nonz, n_nonz_max, nberr)
                   RET_ON_NBERR(nberr)


               enddo
            enddo
         enddo
      enddo
   enddo



   ! The full coupling matrix is dimension n_dof x n_dof
   ! At most n_nonz_max of these can be n_nonzero based on the pairs lying on adjacent elements
   ! We now have a valid CSR indexing of this subset of n_nonz_max pairs

   ! But more of these can be eliminated. Not quite sure why.

   ! Because most meshes will have interesting structure and we conservatively guessed the number
   ! of possible nonzero entries, will likely be zero entries we can remove

   do i=1,this%n_dof-1
      row_lo = this%v_col_ptr(i)
      row_hi = this%v_col_ptr(i+1) - 1

      do j=row_lo,row_hi
         if(row_ind_tmp(j) .ne. 0) cycle   ! this elt is busy, check the next

         ! we've found a zero element in this col all the rest of this col will be zeros, so remove them all

         row_lo2 = this%v_col_ptr(i) + j - row_lo
         this%v_col_ptr(i+1) = row_lo2         ! bring the start of the next col forward
         row_hi2 = this%v_col_ptr(i+2) - 1       ! find the end of that next col
         do k=row_hi+1,row_hi2
            k1 = row_lo2 + k - (row_hi+1)  ! shuffle the columns in that col forward into the empty space
            row_ind_tmp(k1) = row_ind_tmp(k)
            row_ind_tmp(k) = 0
         enddo

         exit  ! this row is done, go back to the next

      enddo


   enddo



   ! squeeze the last col
   i = this%n_dof
   row_lo = this%v_col_ptr(i)
   row_hi = this%v_col_ptr(i+1) - 1
   do j=row_lo,row_hi
      if(row_ind_tmp(j) .eq. 0) then
         row_lo2 = this%v_col_ptr(i) + j - row_lo
         this%v_col_ptr(i+1) = row_lo2
         exit
      endif
   enddo


   ! Now we know n_nonz
   call integer_alloc_1d(this%v_row_ind, n_nonz, 'this%v_row_ind', nberr); RET_ON_NBERR(nberr)

      this%v_row_ind = row_ind_tmp
   this%n_nonz = n_nonz

   ! Find the longest row
   max_col_len = 0
   do i=1,this%n_dof
      row_lo = this%v_col_ptr(i)
      row_hi = this%v_col_ptr(i+1) - 1
      row_len = row_hi - row_lo + 1
      if (row_len .gt. max_col_len) max_col_len = row_len
   enddo

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
