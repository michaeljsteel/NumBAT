#include "numbat_decl.h"

 !  dof_props = 0  => interior ddl (ddl = Degree Of Freedom)
 !  dof_props != 0 => boundary ddl

 !  bdy_cdn = 0 => Dirichlet boundary condition (E-fi_eld: electric wall condition)
 !  bdy_cdn = 1 => Neumann boundary condition (E-fi_eld: magnetic wall condition)
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

subroutine SparseCSC_AC_set_boundary_conditions(this, bdy_cdn, mesh_raw, nberr)

   use numbatmod
   use class_MeshRaw

   class(SparseCSC_AC) :: this
   integer(8) bdy_cdn
   type(MeshRawAC) mesh_raw
   type(NBError) nberr

   ! locals

   integer(8) i_msh, n_dof
   logical is_bdy_mshpt
   integer(8) :: vec3(3)
   vec3 = (/ 1_8, 2_8, 3_8/)


   call integer_alloc_2d(this%m_eqs, 3_8, mesh_raw%n_msh_pts, 'm_eqs', nberr); RET_ON_NBERR(nberr)

   n_dof = 0
   if(bdy_cdn .eq. BCS_DIRICHLET) then   !  all interior points have a degree of freedom

      do i_msh=1,mesh_raw%n_msh_pts
         !is_interior = mesh_raw%v_mshpt_physindex(i) == 0
         is_bdy_mshpt = mesh_raw%is_boundary_mesh_point(i_msh)

         if (.not. is_bdy_mshpt ) then !  each element is associated to 3 interior DOF
            ! Uncommented code implements this:
            !this%m_eqs(1,i_msh) = n_dof + 1
            !this%m_eqs(2,i_msh) = n_dof + 2
            !this%m_eqs(3,i_msh) = n_dof + 3

            this%m_eqs(:,i_msh) = n_dof + vec3
            n_dof = n_dof + 3
         else
            ! Uncommented code implements this:
            !    this%m_eqs(1,i_msh) = 0
            !    this%m_eqs(2,i_msh) = 0
            !    this%m_eqs(3,i_msh) = 0

            this%m_eqs(:,i_msh) = 0
         endif
      enddo

   elseif(bdy_cdn .eq. BCS_NEUMANN) then !  all points have a degree of freedom

      do i_msh=1,mesh_raw%n_msh_pts
         ! Uncommented code implements this:
         ! this%m_eqs(1,i_msh) = n_dof + 1
         ! this%m_eqs(2,i_msh) = n_dof + 2
         ! this%m_eqs(3,i_msh) = n_dof + 3

         this%m_eqs(:,i_msh) = n_dof + vec3
         n_dof = n_dof + 3
      enddo

   endif

   this%n_dof = n_dof
end



subroutine SparseCSC_AC_make_csc_arrays(this, bdy_cdn, mesh_raw, nberr)

   class(SparseCSC_AC) :: this
   type(MeshRawAC) :: mesh_raw

   type(NBError) nberr

   integer(8) bdy_cdn

   ! ------------------------------------------

   integer(8) n_nonz_max, max_col_len, n_nonz

   call this%set_boundary_conditions(bdy_cdn, mesh_raw, nberr); RET_ON_NBERR(nberr)

   this%n_nonz=0

   !Non-sparse matrix is n_dof x n_dof, so  v_col_ptr has length n_dof+1
   call integer_alloc_1d(this%v_col_ptr, this%n_dof+1, 'v_col_ptr', nberr); RET_ON_NBERR(nberr)

   ! sets up a v_col_ptr assuming that every column has the maximum possible non-zero entries
   call this%make_col_ptr_provisional (mesh_raw, n_nonz_max)

   ! v_col_ptr now has the right length for CSC and n_nonz_max is an upper bound for the number of n_nonzeros.
   ! Now get the row_indexes.
   call this%make_arrays_final (mesh_raw, n_nonz_max, n_nonz, max_col_len, nberr);
   RET_ON_NBERR(nberr)


   ! At this point, the row_indices for a given column are in random order
   ! Now we sort them column by column so the rows appear in ascending order within each column
   ! This is another reverse passing to a CSR routine
   call sort_csr (this%n_dof, this%n_nonz, max_col_len, this%v_row_ind, this%v_col_ptr)


   !call this%dump_csc_arrays()

   ! Now we know how big the data arrays are
   call complex_alloc_1d(this%mOp_stiff, this%n_nonz, 'cscmat%mOp_stiff', nberr); RET_ON_NBERR(nberr)

   call complex_alloc_1d(this%mOp_mass, this%n_nonz, 'cscmat%mOp_mass', nberr); RET_ON_NBERR(nberr)


   !  convert from 1-based to 0-based
   !  ----------------------------------------------------------------
   !  Our CSC indexing so far, i.e., ip_col_ptr, has been 1-based
   !  But external calls to functions like valpr.f will need 0-based indexing.
   !  So the indices need shifting by one.
   !  (But not the data arrays: c[0] and fortran[1] refer to the same location, so that just works)

   this%v_row_ind = this%v_row_ind - 1
   this%v_col_ptr = this%v_col_ptr - 1

   return

end subroutine



! subroutine SparseCSC_AC_cscmat_contains_elt_row_col(this, row, col, found, val)

!    class(SparseCSC_AC) :: this
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


subroutine SparseCSC_AC_dump_csc_arrays(this)
   class(SparseCSC_AC) :: this
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

! Find upper bound for number of nonzero elements of FEM operators
! Builds a first version of v_col_ptr
! with enough spaces in each column to allow interactions with every dof on
! neighbouring elements
subroutine SparseCSC_AC_make_col_ptr_provisional (this, mesh_raw, nonz_max)
   use numbatmod

   class(SparseCSC_AC) :: this
   type(MeshRawAC) mesh_raw

   integer(8), intent(out) :: nonz_max

   integer(8) :: i_dof,  i_el, dof, i_nd_xy, tag_mshpt, nd_i
   integer(8) :: last_ct, this_ct

   !TODO: try and use pointer alias for readability
   ! can't make compiler work
   !integer(8), pointer, dimension(:) :: vcp
   !vcp => this%v_col_ptr

   this%v_col_ptr = 0

   !  Determination of the bandwidths


   ! Each dof can play multiple roles according to how many elements it participates in
   ! (eg corner>=3, edge=1 or 2)
   ! Counting these gives an upper bound to the nonzero elements of the FEM matrices

   do i_el=1,mesh_raw%n_msh_elts                          ! for every element
      do nd_i=1,P2_NODES_PER_EL                           ! and all nodes in the element
         tag_mshpt = mesh_raw%m_elnd_to_mshpt(nd_i,i_el)  ! find global mesh point label
         do i_nd_xy=1,3                                    ! for each coordinate
            dof = this%m_eqs(i_nd_xy, tag_mshpt)           ! find index of absolute dof
            if (dof .ne. 0) this%v_col_ptr(dof) = this%v_col_ptr(dof)+1           ! increment number of roles this dof plays in multiple elements
         enddo
      enddo
   enddo

   ! first estimate of nonzeros
   ! for each dof which participates in this%v_col_ptr(i) elts, count every member of 1st, and then one less of the remainder? Explain
   ! each msh_pt is identified as belonging to this%v_col_ptr(dof) elements
   ! each contributes 3*P2_NODES_PER_EL but we avoid double counting _this_ dof by subtracting it's own xyz comps (this%v_col_ptr(i)-1) times
   ! incr = this%v_col_ptr(dof) * 3*P2_NODES_PER_EL - (this%v_col_ptr(dof)-1)*3
   !      = (this%v_col_ptr(dof)-1 +1 ) * 3*P2_NODES_PER_EL - (this%v_col_ptr(dof)-1)*3
   !      = 3*P2_NODES_PER_EL + (this%v_col_ptr(dof)-1 ) * (3*P2_NODES_PER_EL-3)
   !      = 3*P2_NODES_PER_EL + (this%v_col_ptr(dof)-1 ) * 3*(P2_NODES_PER_EL-1)

   nonz_max = 0
   do i_dof=1,this%n_dof
      nonz_max = nonz_max + 3*P2_NODES_PER_EL + 3*(P2_NODES_PER_EL-1)*(this%v_col_ptr(i_dof)-1)
   enddo


   !  Compressed Column Storage (CSC): determine the column pointer

   last_ct = this%v_col_ptr(1)
   this%v_col_ptr(1) = 1
   do i_dof=2,this%n_dof+1
      this_ct = this%v_col_ptr(i_dof)
      this%v_col_ptr(i_dof) = this%v_col_ptr(i_dof-1) + 3*P2_NODES_PER_EL + 3*(P2_NODES_PER_EL-1)*(last_ct-1)
      last_ct = this_ct
   enddo

   nonz_max = this%v_col_ptr(this%n_dof+1) - 1

end



 ! This one is written in CSR format
 ! row/col names seem backward
 ! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine SparseSC_make_arrays_final (this, mesh_raw, n_nonz_max, n_nonz, max_col_len, nberr)

   use numbatmod
   use alloc

   use class_MeshRaw

   class(SparseCSC_AC) :: this
   type(MeshRawAC) mesh_raw

   type(NBError) nberr

   integer(8) n_nonz_max,n_nonz

   integer(8) max_col_len


   integer(8), dimension(:), allocatable :: row_ind_tmp

   integer(8), parameter :: N_ENTITY_PER_EL_AC = 6

   character(len=EMSG_LENGTH) :: emsg
   integer(8) i, j, i_nd, j_nd, k, k1, i_loc_dof, j_loc_dof
   integer(8) i_el, i_dof, i_mshpt, j_dof, j_mshpt
   integer(8) row_lo, row_hi, row_len
   integer(8) row_lo2, row_hi2
   integer(8) ui
   logical found


   ui = stdout

   call integer_alloc_1d(row_ind_tmp, n_nonz_max, 'row_ind_tmp', nberr); RET_ON_NBERR(nberr)

   row_ind_tmp = 0

   ! This code was converted from one in CSR format
   !  Determination of the row indices

   ! For two dof to interact and require nonzero element in the FEM matrices,
   ! they must both exist on the same element

   n_nonz = 0
   do i_el=1,mesh_raw%n_msh_elts

      do i_nd=1,N_ENTITY_PER_EL_AC
         i_mshpt = mesh_raw%m_elnd_to_mshpt(i_nd,i_el)

         do i_loc_dof=1,3
            i_dof = this%m_eqs(i_loc_dof,i_mshpt)

            if (i_dof .eq. 0) cycle  ! not a genuine dof

            row_lo = this%v_col_ptr(i_dof)
            row_hi = this%v_col_ptr(i_dof+1) - 1

            do j_nd=1,N_ENTITY_PER_EL_AC
               j_mshpt = mesh_raw%m_elnd_to_mshpt(j_nd,i_el)

               do j_loc_dof=1,3
                  j_dof = this%m_eqs(j_loc_dof,j_mshpt)
                  if (j_dof .eq. 0) cycle ! not a genuine dof

                  ! Now we know (i_dof, j_dof) can interact and should have a slot
                  ! in the column for i_dof

                  ! Search if the entry (i_dof,j_dof) is already stored
                  found = .false.
                  do k=row_lo,row_hi

                     ! We've run out of nonzero entries in the row indices, and haven't found this element yet
                     ! So create one
                     if(row_ind_tmp(k) .eq. 0) then
                        found = .true.
                        n_nonz =n_nonz + 1

                        if (n_nonz .le. n_nonz_max) then ! We have an empty slot. Go!
                           row_ind_tmp(k) = j_dof
                           exit
                        else  ! should never happen if make_provisional is correct
                           write(emsg, *) "csr_length_AC:n_nonz > n_nonz_max: ",&
                              n_nonz .gt. n_nonz_max
                           call nberr%set(NBERR_BAD_SPARSE_FINAL, emsg)
                           return

                        endif


                     endif

                     ! Entry already exists, nothing to do
                     if(row_ind_tmp(k) .eq. j_dof) then
                        found = .true.
                        exit
                     endif

                  enddo

                  if (found .neqv. .true.) then
                     write(emsg, *) "csr length looking for missing spot"
                     call nberr%set(NBERR_BAD_SPARSE_FINAL, emsg)
                     return
                  endif

               enddo
            enddo
         enddo
      enddo
   end do

   ! Because most meshes will have interesting structure and we conservatively guessed the number
   ! of possible nonzero entries, will likely be zero entries we can remove

   if (n_nonz .lt. n_nonz_max) then ! there are unused entries in some columns

      do i=1,this%n_dof-1
         ! find the initially allocated row slots for this column
         row_lo = this%v_col_ptr(i)
         row_hi = this%v_col_ptr(i+1) - 1

         do j=row_lo,row_hi  ! search for any empties
            if(row_ind_tmp(j) .eq. 0) then ! row indices j to row_hi are empty

               row_lo2 = this%v_col_ptr(i) + j - row_lo ! index of the first nonempty
               this%v_col_ptr(i+1) = row_lo2            ! will be the new first slot in the next column
               row_hi2 = this%v_col_ptr(i+2) - 1        ! current last slot in the next column

               do k=row_hi+1,row_hi2             ! whole next column is pulled back to start at row_lo2
                  k1 = row_lo2 + k - (row_hi+1)
                  row_ind_tmp(k1) = row_ind_tmp(k)
                  row_ind_tmp(k) = 0             ! and are filled in with zeros
               enddo
               !goto 40
               exit
            endif
         enddo
!40       continue
      enddo

      ! now clean up the final column

      i = this%n_dof
      row_lo = this%v_col_ptr(i)
      row_hi = this%v_col_ptr(i+1) - 1

      do j=row_lo,row_hi
         if(row_ind_tmp(j) .eq. 0) then
            row_lo2 = this%v_col_ptr(i) + j - row_lo
            this%v_col_ptr(i+1) = row_lo2  ! Just need to pull down the final column pointer to the last nonzero row
            !goto 50
            exit
         endif
      enddo
!50    continue
   endif

   ! Now we know n_nonz and can set the correct length for the row indices array
   call integer_alloc_1d(this%v_row_ind, n_nonz, 'this%v_row_ind', nberr); RET_ON_NBERR(nberr)
   this%v_row_ind = row_ind_tmp
   this%n_nonz = n_nonz


   ! find maximum row length (actually maximum col length)
   max_col_len = 0
   do i=1,this%n_dof
      row_lo = this%v_col_ptr(i)
      row_hi = this%v_col_ptr(i+1) - 1
      row_len = row_hi - row_lo + 1
      if (row_len .gt. max_col_len) max_col_len = row_len
   enddo



end
