

subroutine MeshAC_allocate(this, n_msh_pts, n_msh_elts, n_elt_mats, nberr)

   class(MeshAC) :: this
   integer(8) :: n_msh_elts, n_msh_pts, n_elt_mats
   type(NBError) nberr

   this%n_msh_pts = n_msh_pts
   this%n_msh_elts = n_msh_elts
   this%n_elt_mats = n_elt_mats

   call double_alloc_2d(this%v_mshpt_xy, 2_8, n_msh_pts, 'v_mshpt_xy', nberr); RET_ON_NBERR(nberr)

   call integer_alloc_1d(this%v_mshpt_physindex, n_msh_pts, 'v_mshpt_physindex', nberr); 
RET_ON_NBERR(nberr)

   call integer_alloc_2d(this%m_elnd_to_mshpt, P2_NODES_PER_EL, n_msh_elts, 'm_elnd_to_mshpt', nberr);
   RET_ON_NBERR(nberr)

end subroutine


 ! Seems identical to the EM version.
subroutine MeshAC_load_mesh_tables_from_py(this, &
   v_mshpt_xy, v_mshpt_physindex, v_elt_material, m_elnd_to_mshpt)

   class(MeshAC) :: this

   integer(8) :: v_mshpt_physindex(this%n_msh_pts)
   double precision :: v_mshpt_xy(2, this%n_msh_pts)
   integer(8) :: v_elt_material(this%n_msh_elts)
   integer(8) :: m_elnd_to_mshpt(P2_NODES_PER_EL, this%n_msh_elts)

   this%v_mshpt_physindex = v_mshpt_physindex
   this%v_mshpt_xy = v_mshpt_xy
   this%v_elt_material = v_elt_material
   this%m_elnd_to_mshpt = m_elnd_to_mshpt

   ! do i=1,this%n_msh_elts
   !    write(*,*) 'eltmat', i, v_elt_material(i)
   ! end do

end subroutine

 ! Seems identical to the EM version.
 ! current not used
subroutine MeshAC_load_mesh_tables_from_scratch(this, mesh_file, dimscale_in_m,nberr)

   class(MeshAC) :: this

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m

   ! outs
   !integer(8) m_elnd_to_mshpt(P2_NODES_PER_EL, this%n_msh_elts)

   type(NBError) nberr

   ! locals
   character(len=EMSG_LENGTH) :: emsg


   double precision xx(2)

   integer(8) n_elt_mats2, n_msh_pts2, n_msh_elts2
   integer(8) i, j, k
   integer(8) ui


   ui = 24

   write(*,*) 'opening meshfile', mesh_file, '<<'

   !  check the mesh file is consistent with what we expect
   open (unit=ui,file=mesh_file, status='old')
   read(ui,*) n_msh_pts2, n_msh_elts2

   if (this%n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", this%n_msh_pts, n_msh_pts2
      call nberr%set(-101_8, emsg)
      return
   endif

   if (this%n_msh_elts .ne. n_msh_elts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_elts != n_msh_elts2 : ", this%n_msh_elts, n_msh_elts2
      call nberr%set(-102_8, emsg)
      return
   endif


   !  Read coordinates of the FEM mesh points
   do i=1,this%n_msh_pts
      read(ui,*) k, (xx(j),j=1,2), this%v_mshpt_physindex(i)
      this%v_mshpt_xy(:,i) = xx*dimscale_in_m
   enddo

   !  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file
   do i=1,this%n_msh_elts

      read(ui,*) k, (this%m_elnd_to_mshpt(j,i),j=1,P2_NODES_PER_EL), this%v_elt_material(i)

      j = this%v_elt_material(i)
      if(n_elt_mats2 .lt. j) n_elt_mats2 = j

      ! shouldn't happen, delete
      if(j .lt. 0) then
         write(emsg,*) "geometry: type_el(i) < 0 : ", i, this%v_elt_material(i)
         call nberr%set(-9_8, emsg)
         return
      endif
   enddo

   close(ui)

   if(n_elt_mats2 .gt. this%n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, this%n_elt_mats
      call nberr%set(-10_8, emsg)
      return
   endif

end


! boundary nodes have non zero GMsh physindex codes
pure logical function  MeshAC_is_boundary_mesh_point(this, msh_pt) result(res)

class(MeshAC), intent(in) :: this
integer(8), intent(in) :: msh_pt

res = this%v_mshpt_physindex(msh_pt) .ne. 0
end function

subroutine MeshAC_find_nodes_for_elt(this, i_el, el_nds_i, el_nds_xy)

   class(MeshAC) :: this
   integer(8) i_el
   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2, P2_NODES_PER_EL)

   integer(8) nd_i, mesh_pt

   do nd_i=1,P2_NODES_PER_EL                          ! For each of the 6 P2 nodes
      mesh_pt = this%m_elnd_to_mshpt(nd_i,i_el)         !    find the global index of that mesh point
      el_nds_i(nd_i) = mesh_pt                        !       store the mesh point indicex
      el_nds_xy(:,nd_i) = this%v_mshpt_xy(:,mesh_pt)  !       and the corresponding physical coordinates
   enddo

end subroutine




