subroutine MeshRaw_allocate(this, n_msh_pts, n_msh_el, n_elt_mats, &
      errco, emsg)

      class(MeshRaw) :: this
      integer(8) :: n_msh_el, n_msh_pts, n_elt_mats
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      this%n_msh_pts = n_msh_pts
      this%n_msh_el = n_msh_el
      this%n_elt_mats = n_elt_mats

      call integer_alloc_1d(this%el_material, n_msh_el, 'el_material', errco, emsg); RETONERROR(errco)

      call double_alloc_2d(this%xy_nodes, 2_8, n_msh_pts, 'xy_nodes', errco, emsg); RETONERROR(errco)

      call integer_alloc_1d(this%node_phys_i, n_msh_pts, 'node_phys_i', errco, emsg); RETONERROR(errco)

      call integer_alloc_2d(this%table_nod, P2_NODES_PER_EL, n_msh_el, 'table_nod', errco, emsg);
      RETONERROR(errco)

   end subroutine

   ! subroutine MeshRaw_destructor(this)
   !    type(MeshRaw) :: this

   ! end subroutine

   subroutine MeshRaw_fill_python_arrays(this, &
      el_material, node_phys_i, table_nod, xy_nodes)

      class(MeshRaw) :: this

      integer(8), intent(out) :: el_material(:)
      integer(8), intent(out) :: node_phys_i(:)
      integer(8), intent(out) :: table_nod(:, :)
      double precision, intent(out) :: xy_nodes(:,:)

      el_material = this%el_material
      node_phys_i = this%node_phys_i
      table_nod = this%table_nod
      xy_nodes = this%xy_nodes

   end subroutine

   ! boundary nodes have non zero GMsh physindex codes
   pure logical function  MeshRaw_is_boundary_node(this, nd) result(res)
      class(MeshRaw), intent(in) :: this
      integer(8), intent(in)  :: nd

      res = this%node_phys_i(nd) .ne. 0
   end function

   pure logical function MeshRaw_is_boundary_node_2(this, i_nd, i_el) result(res)
      class(MeshRaw), intent(in) :: this
      integer(8), intent(in)  :: i_nd, i_el
      res = this%node_phys_i(this%table_nod(i_nd, i_el)) .ne. 0
   end function

   ! get node type by indirection through node table
   integer(8) function  MeshRaw_node_phys_index_by_ref(this, i_nd, i_el) result(res)
      class(MeshRaw) :: this
      integer(8) :: i_nd, i_el
      res = this%node_phys_i(this%table_nod(i_nd, i_el))
   end function






!  Construct the FEM mesh

!  node_phys_i = 0  => interior point
!  node_phys_i != 0 => boundary point specifying which physical line, physindex

!  Reads .mail file to find
!  - x,y coords of mesh points  (xy_nodes)
!  - mesh points associated with each element (table_nod)
!  - whether number of material types read matches expected value n_elt_mats

!  -  Fills:  xy_nodes, node_phys_i, el_material, table_nod

subroutine MeshRaw_construct_node_tables(this, mesh_file, dimscale_in_m, errco, emsg)

   class(MeshRaw) :: this

   ! ins
   character(len=*) mesh_file
   double precision dimscale_in_m

   ! outs


   integer(8) errco
   character(len=EMSG_LENGTH) :: emsg

   ! locals
   double precision xx(2)

   integer(8) n_elt_mats2
   integer(8) n_msh_pts2, n_msh_el2
   integer(8) i, j, k

   integer(8) ui



   ui = 24

   !  check the mesh file is consistent with what we expect
   open (unit = ui, file=mesh_file, status='old')
   read(ui,*) n_msh_pts2, n_msh_el2

   if (this%n_msh_pts .ne. n_msh_pts2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_pts != n_msh_pts2 : ", this%n_msh_pts, n_msh_pts2
      errco=-101
      return
   endif

   if (this%n_msh_el .ne. n_msh_el2) then
      write(emsg,*) "construct_fem_nodal_tables: n_msh_el != n_msh_el2 : ", this%n_msh_el, n_msh_el2
      errco=-102
      return
   endif


!  Read coordinates of the FEM mesh points
   do i=1,this%n_msh_pts
      read(ui,*) k, (xx(j),j=1,2), this%node_phys_i(i)
      this%xy_nodes(:,i) = xx *dimscale_in_m
   enddo

!  Connectivity table
   n_elt_mats2 = 1   !  largest index of materials in the file

   do i=1,this%n_msh_el

      read(ui,*) k, (this%table_nod(j,i),j=1,P2_NODES_PER_EL), this%el_material(i)

      j = this%el_material(i)

      if (n_elt_mats2 .lt. j) n_elt_mats2 = j

   enddo

   close(ui)

   if(n_elt_mats2 .gt. this%n_elt_mats) then
      write(emsg,*) "geometry: n_elt_mats2 > n_elt_mats : ", n_elt_mats2, this%n_elt_mats
      errco=-10
      return
   endif

end



