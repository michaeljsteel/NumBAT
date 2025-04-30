

!#include "numbat_decl.h"


! subroutine SparseCSC_AC_set_boundary_conditions(this, bdy_cdn, mesh_raw,  entities, pbcs, nberr)


!    class(SparseCSC_AC) :: this
!    type(PeriodicBCs) :: pbcs

!    integer(8) :: bdy_cdn !, n_dof
!    integer(8) :: debug

!    type(MeshRaw) :: mesh_raw
!    type(MeshEntities) :: entities

!    type(NBError) nberr


!    !locals
!    double precision, dimension(2,2) :: lat_vecs

!    if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

!       call this%bound_cond_em (bdy_cdn, entities, nberr)

!    elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)
!       debug=0
!       call periodic_lattice_vec (mesh_raw%n_msh_pts, mesh_raw%v_nd_xy, lat_vecs, debug)

!       call periodic_node(mesh_raw%n_msh_el, mesh_raw%n_msh_pts, &
!          P2_NODES_PER_EL, mesh_raw%v_nd_physindex, mesh_raw%v_nd_xy, pbcs%iperiod_N, &
!          pbcs%inperiod_N, mesh_raw%elnd_to_mshpt, lat_vecs)

!       call periodic_N_E_F (entities%n_entities, entities%v_ety_props, entities%v_xy, pbcs%iperiod_N_E_F, &
!          pbcs%inperiod_N_E_F, lat_vecs)

!       call periodic_cond ( bdy_cdn, entities%n_entities, this%n_dof, entities%v_ety_props, &
!          pbcs%iperiod_N_E_F, this%m_eqs, debug)

!    endif

! end subroutine



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

subroutine SparseCSC_AC_set_boundary_conditions(this, bdy_cdn, mesh_raw, nberr)

   use numbatmod
   use class_MeshRaw

   class(SparseCSC_AC) :: this
   integer(8) bdy_cdn
   type(MeshRawAC) mesh_raw
   type(NBError) nberr

   ! locals

   integer(8) i, n_dof
   logical is_interior

   call integer_nalloc_2d(this%m_eqs, 3_8, mesh_raw%n_msh_pts, 'm_eqs', nberr); RET_ON_NBERR(nberr)

   n_dof = 0
   if(bdy_cdn .eq. BCS_DIRICHLET) then   !  all interior points have a degree of freedom

      do i=1,mesh_raw%n_msh_pts
         is_interior = mesh_raw%v_nd_physindex(i) == 0

         if (is_interior ) then !  each element is associated to 3 interior DOF
            this%m_eqs(1,i) = n_dof + 1
            this%m_eqs(2,i) = n_dof + 2
            this%m_eqs(3,i) = n_dof + 3
            n_dof = n_dof + 3
         else
            this%m_eqs(1,i) = 0
            this%m_eqs(2,i) = 0
            this%m_eqs(3,i) = 0
         endif
      enddo

   elseif(bdy_cdn .eq. BCS_NEUMANN) then !  all points have a degree of freedom

      do i=1,mesh_raw%n_msh_pts
         this%m_eqs(1,i) = n_dof + 1
         this%m_eqs(2,i) = n_dof + 2
         this%m_eqs(3,i) = n_dof + 3
         n_dof = n_dof + 3
      enddo

   endif

   this%n_dof = n_dof
end



subroutine SparseCSC_AC_make_csc_arrays(this, mesh_raw, entities, n_nonz_max, n_nonz,  max_row_len, nberr)

   class(SparseCSC_AC) :: this
   type(MeshRawAC) :: mesh_raw
   type(MeshEntitiesAC) :: entities

   type(NBError) nberr

   ! ------------------------------------------

   integer(8) n_nonz_max, max_row_len, n_nonz
   integer(8), dimension(:), allocatable  :: iwork

   !integer(8) row, col, rc_exists, cr_exists, val


   this%n_nonz=0

   call integer_nalloc_1d(this%v_col_ptr, this%n_dof+1, 'v_col_ptr', nberr); RET_ON_NBERR(nberr)


   call this%make_col_ptr_provisional (mesh_raw, n_nonz_max)


   ! v_col_ptr now has the right length for CSC and n_nonz_max is an upper bound for the number of n_nonzeros.
   ! Now get the row_indexes.
   call this%make_arrays_final (mesh_raw, n_nonz_max, n_nonz, max_row_len, nberr);
   RET_ON_NBERR(nberr)






   ! At this point, the row_indices for a given column are in random order
   ! Now we sort them column by column so the rows appear in ascending order within each column
   ! This is another reverse passing to a CSR routine
   call integer_nalloc_1d(iwork, 3*mesh_raw%n_msh_pts, 'iwork', nberr); RET_ON_NBERR(nberr)

   call sort_csr (this%n_dof, this%n_nonz, max_row_len, this%v_row_ind, this%v_col_ptr,  iwork)


   !call this%dump_csc_arrays()

   ! Now we know how big the data arrays are
   call complex_nalloc_1d(this%mOp_stiff, this%n_nonz, 'cscmat%mOp_stiff', nberr); RET_ON_NBERR(nberr)

   call complex_nalloc_1d(this%mOp_mass, this%n_nonz, 'cscmat%mOp_mass', nberr); RET_ON_NBERR(nberr)


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



subroutine SparseCSC_AC_cscmat_contains_elt_row_col(this, row, col, found, val)

   class(SparseCSC_AC) :: this
   integer(8) row, col
   integer(8) found

   integer(8) val
   integer(8) vct_1, vct_2, ri

   found = 0
   if (col .lt. 1 .or. col .gt. this%n_dof) then
      return
   endif

   vct_1 = this%v_col_ptr(col)
   vct_2 = this%v_col_ptr(col+1)
   do ri = vct_1, vct_2
      if (this%v_row_ind(ri) .eq. row) then
         found = 1
         val = 1
      endif
   enddo

end subroutine


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



subroutine SparseCSC_AC_make_col_ptr_provisional (this, mesh_raw,  nonz_max)
   class(SparseCSC_AC) :: this

   type(MeshRawAC) mesh_raw


   integer(8) :: lb(this%n_dof+1)
   integer(8), intent(out) :: nonz_max


   integer(8), parameter :: nddl_0 = 6 !  Different to em case. Different FEM?

   integer(8) :: i, k, iel, ind_ip, ip
   integer(8) :: k_copy1, k_copy2

  lb = 0


   !  Determination of the bandwidths


   do iel=1,mesh_raw%n_msh_el
      do i=1,nddl_0
         ip = mesh_raw%elnd_to_mshpt(i,iel)
         do k=1,3
            ind_ip = this%m_eqs(k,ip)
            if (ind_ip .ne. 0) lb(ind_ip) = lb(ind_ip)+1
         enddo
      enddo
   enddo

   nonz_max = 0
   do i=1,this%n_dof
      nonz_max = nonz_max + 3*nddl_0 + 3*(nddl_0-1)*(lb(i)-1)
   enddo


   !  Compressed Row Storage (CRS): determine the row pointer

   k_copy1 = lb(1)
   lb(1) = 1
   do i=2,this%n_dof+1
      k_copy2 = lb(i)
      lb(i) = lb(i-1) + 3*nddl_0 + 3*(nddl_0-1)*(k_copy1-1)
      k_copy1 = k_copy2
   enddo

   nonz_max = lb(this%n_dof+1) - 1

   this%v_col_ptr=lb

   return
end




 ! subroutine SparseSC_make_arrays_final (this, mesh_raw, entities, n_nonz_max, max_row_len, errco, emsg)


 !    class(SparseCSC_AC) :: this
 !    type(MeshRaw) :: mesh_raw
 !    type(MeshEntities) :: entities

 !    integer(8) n_nonz_max, max_row_len
 !    integer(8) errco
 !    character(len=EMSG_LENGTH) emsg

 !    ! csr_length labels v_row_ind and v_col_ptr in reverse to here!
 !    ! length of v_row_ind is determined inside csr_length and so allocated there

 !    write(*,*) 'maf 1'

 !    call csr_length (mesh_raw%n_msh_el, entities%n_entities, this%n_dof,  entities%v_tags, this%m_eqs, &
 !       this%v_row_ind, this%v_col_ptr, n_nonz_max, this%n_nonz, max_row_len, errco, emsg)
 !    RETONERROR(errco)



 !end subroutine



 ! This one is written in CSR format
 ! row/col names seem backward
 ! this seems to be a row-like csr converted to a column-like csr with no name changes?

subroutine SparseSC_make_arrays_final (this, mesh_raw, n_nonz_max, n_nonz, max_row_len, nberr)

   use numbatmod
   use alloc

   use class_MeshRaw

   class(SparseCSC_AC) :: this
   type(MeshRawAC) mesh_raw

   type(NBError) nberr

   integer(8) n_nonz_max,n_nonz




   integer(8) max_row_len


   integer(8), dimension(:), allocatable :: row_ind_tmp

   integer(8), parameter :: N_ENTITY_PER_EL_AC = 6

   integer(8) i, j, i_nd, j_nd, k, k1, i_dof, j_dof
   integer(8) iel, ind_ip, ip, ind_jp, jp
   integer(8) row_start, row_end, row_len
   integer(8) row_start2, row_end2
   integer(8) ui


   ui = stdout

   call integer_nalloc_1d(row_ind_tmp, n_nonz_max, 'row_ind_tmp', nberr); RET_ON_NBERR(nberr)

   row_ind_tmp = 0

   ! This code was converted from one in CSR format
   !  Determination of the row indices

   n_nonz = 0
   do iel=1,mesh_raw%n_msh_el

      do i_nd=1,N_ENTITY_PER_EL_AC
         ip = mesh_raw%elnd_to_mshpt(i_nd,iel)

         do i_dof=1,3
            ind_ip = this%m_eqs(i_dof,ip)
            if (ind_ip .ne. 0) then
               row_start = this%v_col_ptr(ind_ip)
               row_end = this%v_col_ptr(ind_ip+1) - 1

               do j_nd=1,N_ENTITY_PER_EL_AC
                  jp = mesh_raw%elnd_to_mshpt(j_nd,iel)

                  do j_dof=1,3
                     ind_jp = this%m_eqs(j_dof,jp)
                     if (ind_jp .ne. 0) then
                        ! Search if the entry (ind_ip,ind_jp) is already stored
                        do k=row_start,row_end
                           if(row_ind_tmp(k) .eq. 0) goto 20
                           if(row_ind_tmp(k) .eq. ind_jp) goto 30
                        enddo

                        print*, "csr_length_AC: There is a problem!", " Aborting..."
                        stop

                        20                      continue

                        !  No entry exists for (ind_ip,ind_jp); create new one
                        n_nonz =n_nonz + 1
                        if (n_nonz .gt. n_nonz_max) then
                           print*, "csr_length_AC:n_nonz > n_nonz_max: ",&
                           n_nonz .gt. n_nonz_max
                           print*, "csr_length_AC: Aborting..."
                           stop
                        endif

                        row_ind_tmp(k) = ind_jp
                        30                      continue

                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
   end do



   ! squeeze away the zero entries
   ! added so as to handle more type of domains/meshes

   if (n_nonz .lt. n_nonz_max) then
      do i=1,this%n_dof-1
         row_start = this%v_col_ptr(i)
         row_end = this%v_col_ptr(i+1) - 1
         do j=row_start,row_end
            if(row_ind_tmp(j) .eq. 0) then
               row_start2 = this%v_col_ptr(i) + j - row_start
               this%v_col_ptr(i+1) = row_start2
               row_end2 = this%v_col_ptr(i+2) - 1
               do k=row_end+1,row_end2
                  k1 = row_start2 + k - (row_end+1)
                  row_ind_tmp(k1) = row_ind_tmp(k)
                  row_ind_tmp(k) = 0
               enddo
               goto 40
            endif
         enddo
40       continue
      enddo
      i = this%n_dof
      row_start = this%v_col_ptr(i)
      row_end = this%v_col_ptr(i+1) - 1
      do j=row_start,row_end
         if(row_ind_tmp(j) .eq. 0) then
            row_start2 = this%v_col_ptr(i) + j - row_start
            this%v_col_ptr(i+1) = row_start2
            goto 50
         endif
      enddo
50    continue
   endif

   max_row_len = 0
   do i=1,this%n_dof
      row_start = this%v_col_ptr(i)
      row_end = this%v_col_ptr(i+1) - 1
      row_len = row_end - row_start + 1
      if (row_len .gt. max_row_len) max_row_len = row_len
   enddo

   write(*,*) 'final', max_row_len

   ! if ((ipointer+n_nonz) .gt. int_max) then
   !    write(ui,*) "csr_length_AC: (ipointer+nonz) > int_max : ",&
   !    &(ipointer+n_nonz), int_max
   !    write(ui,*) "csr_length_AC: n_nonz_max = ", n_nonz_max
   !    write(ui,*) "csr_length_AC: increase the size of int_max"
   !    write(ui,*) "csr_length_AC: Aborting..."
   !    stop
   ! else



   ! Now we know n_nonz
   call integer_nalloc_1d(this%v_row_ind, n_nonz, 'this%v_row_ind', nberr); RET_ON_NBERR(nberr)

   this%v_row_ind = row_ind_tmp


      ! !  Copy the local array row_ind_tmp into col_ind
      ! do i=1,n_nonz
      !    col_ind(i) = row_ind_tmp(i)
      ! enddo
   ! endif

   this%n_nonz = n_nonz

end
