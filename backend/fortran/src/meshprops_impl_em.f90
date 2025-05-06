subroutine MeshRawEM_allocate(this, n_msh_pts, n_msh_elts, n_elt_mats, nberr)

   class(MeshRawEM) :: this
   integer(8) :: n_msh_elts, n_msh_pts, n_elt_mats
   type(NBError) nberr

   this%n_msh_pts = n_msh_pts
   this%n_msh_elts = n_msh_elts
   this%n_elt_mats = n_elt_mats

   call integer_alloc_1d(this%v_elt_material, n_msh_elts, 'v_elt_material', nberr); RET_ON_NBERR(nberr)

   call double_alloc_2d(this%v_mshpt_xy, 2_8, n_msh_pts, 'v_mshpt_xy', nberr); RET_ON_NBERR(nberr)

   call integer_alloc_1d(this%v_mshpt_physindex, n_msh_pts, 'v_mshpt_physindex', nberr); RET_ON_NBERR(nberr)

   call integer_alloc_2d(this%m_elnd_to_mshpt, P2_NODES_PER_EL, n_msh_elts, 'm_elnd_to_mshpt', nberr);
   RET_ON_NBERR(nberr)

end subroutine

 ! subroutine MeshRawEM_destructor(this)
 !    type(MeshRawEM) :: this

 ! end subroutine

subroutine MeshRawEM_fill_python_arrays(this, &
   v_elt_material, v_mshpt_physindex, m_elnd_to_mshpt, v_mshpt_xy)

   class(MeshRawEM) :: this

   integer(8), intent(out) :: v_elt_material(:)
   integer(8), intent(out) :: v_mshpt_physindex(:)
   integer(8), intent(out) :: m_elnd_to_mshpt(:, :)
   double precision, intent(out) :: v_mshpt_xy(:,:)

   v_elt_material = this%v_elt_material
   v_mshpt_physindex = this%v_mshpt_physindex
   m_elnd_to_mshpt = this%m_elnd_to_mshpt
   v_mshpt_xy = this%v_mshpt_xy

end subroutine

 ! boundary nodes have non zero GMsh physindex codes
pure logical function  MeshRawEM_is_boundary_mesh_point(this, msh_pt) result(res)
   class(MeshRawEM), intent(in) :: this
   integer(8), intent(in) :: msh_pt

   res = this%v_mshpt_physindex(msh_pt) .ne. 0
end function

pure logical function MeshRawEM_is_boundary_node_at_element(this, i_nd, i_el) result(res)
   class(MeshRawEM), intent(in) :: this
   integer(8), intent(in)  :: i_nd, i_el
   res = this%v_mshpt_physindex(this%m_elnd_to_mshpt(i_nd, i_el)) .ne. 0
end function

 ! get node type by indirection through node table
integer(8) function  MeshRawEM_node_phys_index_by_ref(this, i_nd, i_el) result(res)
   class(MeshRawEM) :: this
   integer(8) :: i_nd, i_el
   res = this%v_mshpt_physindex(this%m_elnd_to_mshpt(i_nd, i_el))
end function






!  Construct the FEM mesh

!  v_mshpt_physindex = 0  => interior point
!  v_mshpt_physindex != 0 => boundary point specifying which physical line, physindex

!  Reads .mail file to find
!  - x,y coords of mesh points  (v_mshpt_xy)
!  - mesh points associated with each element (m_elnd_to_mshpt)
!  - whether number of material types read matches expected value n_elt_mats

!  -  Fills:  v_mshpt_xy, v_mshpt_physindex, v_elt_material, m_elnd_to_mshpt

subroutine MeshRawEM_construct_mesh_tables(this, mesh_file, dimscale_in_m, nberr)

   class(MeshRawEM) :: this

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m

   ! outs
   type(NBError) nberr


   ! locals
   character(len=EMSG_LENGTH) :: emsg
   double precision xx(2)

   integer(8) n_elt_mats_max
   integer(8) n_msh_pts2, n_msh_elts2
   integer(8) i, j, k

   integer(8) ui



   ui = 24

   !  check the mesh file is consistent with what we expect
   open (unit = ui, file=mesh_file, status='old')
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
      this%v_mshpt_xy(:,i) = xx *dimscale_in_m
   enddo

!  Connectivity table
   n_elt_mats_max = 1   !  largest index of materials in the file

   do i=1,this%n_msh_elts

      read(ui,*) k, (this%m_elnd_to_mshpt(j,i),j=1,P2_NODES_PER_EL), this%v_elt_material(i)

      j = this%v_elt_material(i)

      if (n_elt_mats_max .lt. j) n_elt_mats_max = j

   enddo

   close(ui)

   if(n_elt_mats_max .gt. this%n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats_max > n_elt_mats : ", n_elt_mats_max, this%n_elt_mats
      call nberr%set(-10_8, emsg)
      return
   endif

end



subroutine MeshRawEM_find_nodes_for_elt(this, i_el, &
   el_nds_i, el_nds_xy, is_curved)

   class(MeshRawEM) :: this
   integer(8) i_el
   integer(8) el_nds_i(P2_NODES_PER_EL)
   double precision el_nds_xy(2,P2_NODES_PER_EL)
   logical is_curved

   integer(8) nd_i, mesh_pt

   do nd_i=1,P2_NODES_PER_EL                    ! For each of the 6 P2 nodes
      mesh_pt = this%m_elnd_to_mshpt(nd_i,i_el)    !    find the index of the mesh point
      el_nds_i(nd_i) = mesh_pt                      !    store the mesh point indices for this element
      el_nds_xy(:,nd_i) = this%v_mshpt_xy(:,mesh_pt)  !    find their physical positions
   enddo

   is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, el_nds_xy)

end subroutine












