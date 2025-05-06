
subroutine csc_make_col_ptr_provisional_AC (mesh_raw, cscmat,   nonz_max)

   use numbatmod

   use class_MeshRaw
   use class_SparseCSC_AC

   type(MeshRawAC) mesh_raw
   type(MeshEntitiesAC) entities
   type(SparseCSC_AC) cscmat


   !integer(8), intent(in) :: n_msh_elts, n_msh_pts
   !integer(8), intent(in) ::  elnd_to_mshpt (P2_NODES_PER_EL,mesh_raw%n_msh_elts)

   integer(8) :: lb(cscmat%n_dof+1)
   integer(8), intent(out) :: nonz_max


   !  Local variables
   integer(8), parameter :: nddl_0 = 6 !  Different to em case. Different FEM?

   integer(8) :: i, k, iel, ind_ip, ip
   integer(8) :: k_copy1, k_copy2



   !do i=1,n_dof+1
   !lb(i) = 0
   !enddo
   lb = 0


   !  Determination of the bandwidths


   do iel=1,mesh_raw%n_msh_elts
      do i=1,nddl_0
         ip = mesh_raw%elnd_to_mshpt(i,iel)
         do k=1,3
            ind_ip = cscmat%m_eqs(k,ip)
            if (ind_ip .ne. 0) lb(ind_ip) = lb(ind_ip)+1
         enddo
      enddo
   enddo

   nonz_max = 0
   do i=1,cscmat%n_dof
      nonz_max = nonz_max + 3*nddl_0 + 3*(nddl_0-1)*(lb(i)-1)
   enddo


   !  Compressed Row Storage (CRS): determine the row pointer


   k_copy1 = lb(1)
   lb(1) = 1
   do i=2,cscmat%n_dof+1
      k_copy2 = lb(i)
      lb(i) = lb(i-1) + 3*nddl_0 + 3*(nddl_0-1)*(k_copy1-1)
      k_copy1 = k_copy2
   enddo

   nonz_max = lb(cscmat%n_dof+1) - 1

   cscmat%v_col_ptr=lb

   return
end
