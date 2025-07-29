
subroutine MeshEntities_allocate(this, n_msh_elts, nberr)

   class(MeshEntities) :: this
   integer(8) :: n_msh_elts
   type(NBError) nberr

   integer(8) max_est_entities

   max_est_entities= 9 * n_msh_elts

   call integer_alloc_2d(this%v_tags, 14_8, n_msh_elts, 'v_tags', nberr); RET_ON_NBERR(nberr)
   call double_alloc_2d(this%v_xy, 2_8, max_est_entities, 'v_xy', nberr); RET_ON_NBERR(nberr)
   call integer_alloc_2d(this%v_ety_props, 2_8, max_est_entities, 'x_ety_props', nberr);
   RET_ON_NBERR(nberr)


   call integer_alloc_1d(this%visited, max_est_entities, 'visited', nberr); RET_ON_NBERR(nberr)

   !  Define endpoints of the 3 edges (mid-point) of the reference triangle

   !i = 1
   this%edge_end_nodes(1,1) = 1
   this%edge_end_nodes(2,1) = 2

   !i = 2
   this%edge_end_nodes(1,2) = 2
   this%edge_end_nodes(2,2) = 3

   !i = 3
   this%edge_end_nodes(1,3) = 1
   this%edge_end_nodes(2,3) = 3


end subroutine


 ! !  Storage locations in sequence
 ! !  - tab_N_E_F = table_N_E_F,   shape: 14 x n_msh_elts
 ! !  - table_edges     shape: 4 x n_msh_pts
 ! !

 ! !  V = number of vertices
 ! !  E = number of edges
 ! !  F = number of faces
 ! !  C = number of cells (3D, triangle)
 ! !
 ! !  From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
 ! !  n_msh_pts = (number of vertices) + (number of mid-edge point) = V + E;
 ! !

subroutine MeshEntities_build_mesh_tables(this, mesh, nberr)

   use alloc

   class(MeshEntities) :: this
   type(MeshEM) :: mesh

   type(NBError) nberr

   ! locals
   integer(8) ui_out

   !  ----------------------------------------------

   ui_out = stdout

   call this%allocate(mesh%n_msh_elts, nberr)
   RET_ON_NBERR(nberr)

   call this%check_bdy_elements_are_consistent(mesh, nberr)
   RET_ON_NBERR(nberr)


   !  Fills:  v_tags[1,:] = [1..n_msh_elts]
   call this%count_faces (mesh%n_msh_elts)

   !  Fills: n_edge, v_tags[2:4,:]
   call this%count_edges (mesh)

   !  Fills: v_tags[5:,:], visited[1:n_msh_pts], n_msh_pts_3
   call this%count_nodes_P3 (mesh, nberr); RET_ON_NBERR(nberr)

   ! Every distinct entity has been uniquely tagged
   ! Total number of labelled structure is now known
   this%n_entities = this%n_edges + this%n_faces + this%n_msh_pts_p3


   !  Fill v_xy and v_ety_props arrays
   call this%analyse_face_and_edges (mesh)
   call this%analyse_p3_nodes(mesh)

   deallocate(this%visited)

end subroutine

! Not sure why this is an issue
! Could be broken .gmsh template file with bad physical line?
subroutine MeshEntities_check_bdy_elements_are_consistent(this, mesh, nberr)

   class(MeshEntities) :: this
   type(MeshEM) :: mesh

   type(NBError) nberr

   character(len=EMSG_LENGTH) :: emsg

   integer(8) el_i, edge_nd
   integer(8) ed_vert_nda, ed_vert_ndb  ! edge vertices nodes
   ! ---------------------------------------------

   do el_i=1,mesh%n_msh_elts  ! for each element

      ! Check that boundary elements are well constructed.
      ! for its edge nodes 4,5,6
      do edge_nd=4, P2_NODES_PER_EL

         ! if edge node is on a physical bdy
         if (mesh%is_boundary_mesh_point_by_elt_node(edge_nd,el_i)) then

            ! find node indices (1,2,3) of vertices of the edge
            ed_vert_nda = this%edge_end_nodes(1, edge_nd-3)
            ed_vert_ndb = this%edge_end_nodes(2, edge_nd-3)

            ! Check that both vertices are also bdy points
            ! (else would be a broken mesh)

            if (.not. mesh%is_boundary_mesh_point_by_elt_node(ed_vert_nda,el_i) &
               .or. .not. mesh%is_boundary_mesh_point_by_elt_node(ed_vert_ndb,el_i)) then

               write(emsg,*) "list_edge: v_tags = ", &
                  mesh%node_phys_index_by_ref(ed_vert_nda,el_i), &
                  mesh%node_phys_index_by_ref(ed_vert_ndb,el_i), &
                  mesh%node_phys_index_by_ref(edge_nd,el_i), &
                  "node_phys_i(ed_vert_nda) = ", mesh%m_elnd_to_mshpt(ed_vert_nda,el_i), &
                  "node_phys_i(ed_vert_ndb) = ", mesh%m_elnd_to_mshpt(ed_vert_ndb,el_i), &
                  "node_phys_i(edge_nd) = ", mesh%m_elnd_to_mshpt(edge_nd,el_i)
               call nberr%set(NBERR_BAD_MESH_EDGES, emsg)
               return

            endif
         endif
      enddo
   enddo

end subroutine


! fills v_tags[1, 1_n_msh_elts]
subroutine MeshEntities_count_faces (this, n_elts)

   class(MeshEntities) :: this
   integer(8) n_elts, i

   ! Every element has one face
   this%n_faces = n_elts

   !  The number of each face is the number of its owning element
   do i=1,n_elts
      this%v_tags(1,i) = i
   enddo

end


! Counts the edges and fills v_tags[2..4, 1_n_msh_elts]
!   corresponding to edges 1,2,3 of each element.
! Each edge is assigned a unique edge label starting from n_msh_elts+1
!
! Where an edge lies between two elements it gets the same label
! based on whichever is found first.
! v_tags[ety_i, el_i] = v_tags[ety_j, el_j] if they are the same mesh point

subroutine MeshEntities_count_edges (this, mesh)
   class(MeshEntities) :: this
   type(MeshEM) :: mesh

   ! -----------------------------------------------
   integer(8) n_edge
   integer(8) el_i, edge_nd, edge_num, edge_ety
   !integer(8) ed_vert_nda, ed_vert_ndb
   integer(8) nd_mshpt
   integer(8) new_lab, old_lab
   ! ---------------------------------------------

   this%visited= 0
   n_edge = 0


   do el_i=1,mesh%n_msh_elts  ! for each element

      do edge_num=1,3

         edge_nd = edge_num+3 ! is nd 4,5,6
         edge_ety = edge_num+1  ! edges 1,2,3 are stored as mesh_entities 2,3,4

         nd_mshpt = mesh%m_elnd_to_mshpt(edge_nd, el_i)  ! find the node

         old_lab = this%visited(nd_mshpt)

         if (old_lab .eq. 0) then        ! a new edge encountered
            n_edge = n_edge + 1

            new_lab = n_edge+ this%n_faces

            ! find node numbers 1,2,or3 of the two ends of this edge
            ! unused?
            !ed_vert_nda = this%edge_end_nodes(1, edge_num)
            !ed_vert_ndb = this%edge_end_nodes(2, edge_num)

            this%v_tags(edge_ety, el_i) = new_lab

            this%visited(nd_mshpt) = new_lab  ! visited stores its edge number

         else  ! Met before, take the
            this%v_tags(edge_ety,el_i) = old_lab
         endif

      enddo
   enddo

   this%n_edges = n_edge

end


! get absolute mesh points at elt el_i for vertex nodes surrounding edge node nd_i=4,5,6
subroutine find_vertex_mesh_points_of_edge(mesh, el_i, nd_i, vertices)
   type(MeshEM) :: mesh
   integer(8) el_i, nd_i
   integer(8) nd_vert
   integer(8) vertices(2)

   nd_vert = nd_i - 3
   vertices(1) = mesh%m_elnd_to_mshpt(nd_vert, el_i)

   nd_vert = nd_i - 2
   if (nd_vert .gt. 3) nd_vert = nd_vert - 3
   vertices(2) = mesh%m_elnd_to_mshpt(nd_vert,el_i)

end subroutine


! v_tags[5..7,:] - labels for P3 vertices
! unique tags keep counting from n_faces + n_edges+1
! 10 P3 nodes in order are 3 vertices, six 1/3 + 2/3 edge nodes, and 1 interior node
subroutine MeshEntities_count_nodes_P3 (this, mesh, nberr)

   class(MeshEntities) :: this
   type(MeshEM) :: mesh


   type(NBError) nberr

   ! Locals
   character(len=EMSG_LENGTH) :: emsg
   integer(8) elt_mshpts(P2_NODES_PER_EL)
   integer(8) nd_j, j1, m1, m2, j
   integer(8) el_i, nd_i, el_i2, nd_i2, nd_mshpt, nd_mshpt2
   integer(8) row_off, lab, lab_off
   integer(8) vert_1(2), vert_2(2)
   integer(8) row_off_vert, row_off_P3_edges, row_off_P3_interior

   this%visited = 0

   ! The first 4 entries of this%v_tags(*,i) correspond to face and edges
   ! These have already been filed in count_faces and count_edges

   row_off = 4                             ! filled rows in v_tags
   row_off_vert = 4
   row_off_P3_edges = 7
   row_off_P3_interior = 13

   lab_off = this%n_edges + this%n_faces  ! assigned labels
   lab = 0

   do el_i=1,mesh%n_msh_elts

      ! find the absolute node indices of this element
      elt_mshpts = mesh%m_elnd_to_mshpt(:, el_i)

      !  P3 element: Vertices
      do nd_i=1,3
         nd_mshpt = elt_mshpts(nd_i)          ! find absolute node index

         if (this%visited(nd_mshpt) .eq. 0) then  ! This vertex is new
            this%visited(nd_mshpt) = el_i         !  tag it with the owning elt

            lab = lab + 1                         ! assign a new tag to this vertex
            this%v_tags(nd_i+row_off,el_i) = lab + lab_off

         else                                     ! find and copy the already assigned tag
            el_i2 = this%visited(nd_mshpt)        ! find the owning elt

            ! find which vertex node of the owning el_i2 elt has a matching mesh point
            nd_i2 = 0
            do nd_j=1,3
               nd_mshpt2 = mesh%m_elnd_to_mshpt(nd_j, el_i2)
               if (nd_mshpt .eq. nd_mshpt2) nd_i2 = nd_j
            enddo

            ! illformed mesh
            if(nd_i2 .eq. 0) then
               write(emsg, *) "list_node_P3: problem with a vertex ", el_i, nd_i, &
                  "nd, visited(nd) = ", nd_mshpt, this%visited(nd_mshpt)
               call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
               return
            endif

            ! copy the vertex label
            this%v_tags(nd_i+row_off,el_i) = this%v_tags(nd_i2+row_off,el_i2)
         endif
      enddo


      !  Have filled 3 more rows, row_off gets +3

      !  P3 element: 6 1/3-2/3 nodes on the edges surrounding the P2 mid-edge nodes
      do nd_i=4,6                       ! for each P2 edge node
         nd_mshpt = elt_mshpts(nd_i)

         !  Find absolute mesh points of the vertices of that edge
         call find_vertex_mesh_points_of_edge(mesh, el_i, nd_i, vert_1)

         if (this%visited(nd_mshpt) .eq. 0) then  ! new edge
            this%visited(nd_mshpt) = el_i         ! tag it with the owning element
            m1 = 2*(nd_i-4)+1                     ! identify its end vertices as 1.2, 3.4, or 5.6
            m2 = m1+1
            do j1=m1,m2
               lab = lab+1
               this%v_tags(j1+row_off_P3_edges,el_i) = lab + lab_off
            enddo

         else
            el_i2 = this%visited(nd_mshpt)        ! find the owning elt
            nd_i2 = 0

            ! find which edge node of the owning elt matches
            do nd_j=4,6
               nd_mshpt2=mesh%m_elnd_to_mshpt(nd_j,el_i2)

               if (nd_mshpt .eq. nd_mshpt2) then ! found it
                  nd_i2 = nd_j
                  call find_vertex_mesh_points_of_edge(mesh, el_i2, nd_i2, vert_2)
               endif
            enddo

            if(nd_i2 .eq. 0) then
               write(emsg, *)  "list_node_P3: problem with a node ", el_i, nd_i
               call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
            endif

            ! copy the labels
            do j=1,2
               !  local numbering along the edges
               if (vert_2(1) .eq. vert_1(1) .and. vert_2(2) .eq. vert_1(2)) then
                  !  The nodes on the edges nd_i and nd_i2 are numbered in the same order
                  !  This is possible only when the elements el_i and el_i2 have opposite orientations
                  m1 = j+2*(nd_i-4)
                  m2 = j+2*(nd_i2-4)
               elseif (vert_2(1) .eq. vert_1(2) .and.  vert_2(2) .eq. vert_1(1)) then
                  !  The nodes on the edges nd_i and nd_i2 are numbered in the opposite order
                  !  the numbering of the nodes are reversed
                  j1 = 3 - j
                  m1 = j1+2*(nd_i-4)
                  m2 = j+2*(nd_i2-4)
               else
                  write(emsg,*) "list_node_P3: problems: ", &
                     "Check the edge endpoints", &
                     "nd_i, m_elnd_to_mshpt(nd_i,el_i) = ", nd_i, &
                     mesh%m_elnd_to_mshpt(nd_i,el_i), &
                     "nd_i2, m_elnd_to_mshpt(nd_i2,el_i2) = ", nd_i2, &
                     mesh%m_elnd_to_mshpt(nd_i2,el_i2), &
                     "el_i, el_i2 = ", el_i, el_i2, &
                     "vert_1 = ", vert_1, &
                     "vert_2 = ", vert_2
                  call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
                  return
               endif

               this%v_tags(m1+row_off_P3_edges,el_i) = this%v_tags(m2+row_off_P3_edges,el_i2)
            enddo

         endif
      enddo


      !  Numbering the interior nodes of the triangle
      !  there is only one interior node for a P3 triangle and no possibility of appearing in two elts
      lab = lab+1
      this%v_tags( 1 + row_off_P3_interior, el_i) = lab + lab_off

   enddo

   this%n_msh_pts_p3 = lab

end


! Note that P2 vertices are not explicitly accounted for.
! They seem to be treated as 'faces' and then we get 3 DOF in sparse_csc_em_impl.f90
subroutine MeshEntities_analyse_face_and_edges (this, mesh)


   class(MeshEntities) :: this
   type(MeshEM) :: mesh


   double precision, parameter :: one_third = 1.d0/3.d0
   integer(8) el_i, j, tag, nd_j
   integer(8) phystype_nds(10)
   double precision el_nds_xy(2,6)
   integer(8) el_nds_i(P2_NODES_PER_EL)
   logical is_curved

   this%v_ety_props = 0
   this%visited= 0


   do el_i=1,mesh%n_msh_elts

      call mesh%find_nodes_for_elt(el_i, el_nds_i, el_nds_xy, is_curved)
      do j=1,P2_NODES_PER_EL
         phystype_nds(j) = mesh%v_mshpt_physindex(el_nds_i(j))
      enddo

      ! Face properties
      tag = this%v_tags(ETY_TAG_OFFSET_FACE + 1, el_i)

      ! Position is the barycentre
      this%v_xy(:,tag) = (el_nds_xy(:,1) + el_nds_xy(:,2) + el_nds_xy(:,3)) * one_third

      this%v_ety_props(ETY_PROP_PHYSTYPE, tag) = 0   ! A face is considered interior
      this%v_ety_props(ETY_PROP_DIMENSION, tag) = 2  ! Faces are dimension 2


      !  P2 edge properties
      do j=1,3
         nd_j = j+3   ! P2 eges are nodes 4,5,6
         tag = this%v_tags(ETY_TAG_OFFSET_P2_EDGES + j, el_i)  ! edges are at OP2EGES + 1,2,3


         if (this%visited(tag) .eq. 0) then  ! only do each tag once
            this%visited(tag) = 1

            this%v_xy(:,tag) = el_nds_xy(:,nd_j)

            this%v_ety_props(ETY_PROP_PHYSTYPE, tag) = phystype_nds(nd_j)
            this%v_ety_props(ETY_PROP_DIMENSION, tag) = 1     !  edges are dimension 1
         endif

      enddo
   enddo
end




subroutine MeshEntities_analyse_p3_nodes(this, mesh)

   class(MeshEntities) :: this
   type(MeshEM) :: mesh

   integer(8)  k1, ip(2,3), tag
   integer(8) el_i, nd_i,  nd_i2, ety_off, nd_tag
   integer(8) el_node_tags(P2_NODES_PER_EL)
   integer(8) p3_tags(N_ENTITY_PER_EL)

   double precision xx1, xx2, xx3, yy1, yy2, yy3
   double precision dx1, dy1

   double precision, parameter :: one_third = 1.0d0/3.0d0


!  ip(1,i) = i+1 MOD 3 ; Number of the next vertex to vertex i
!  ip(2,i) = i+2 MOD 3 ; Number of the second next vertex to vertex i
   ip(1,1) = 2
   ip(1,2) = 3
   ip(1,3) = 1

   ip(2,1) = 3
   ip(2,2) = 1
   ip(2,3) = 2


   this%visited = 0


   !  The first 4 entries of this%v_tags(*,i) correspond to face and P2 edges and have been done
   ety_off = 4
   do el_i=1, mesh%n_msh_elts

      el_node_tags = mesh%m_elnd_to_mshpt(:, el_i)

      !  the 10 node of a P3 element
      ! do nd_i=1,P3_NODES_PER_EL
      !    p3_tags(nd_i) = this%v_tags(nd_i+ety_off,el_i)
      ! enddo
      !p3_tags(1:P3_NODES_PER_EL)  = this%v_tags(ety_off+1:ety_off+10, el_i)

      p3_tags(1:P3_NODES_PER_EL)  = this%v_tags(ety_off+1: ety_off+P3_NODES_PER_EL, el_i)



      !  first scan the vertices
      do nd_i=1,3
         nd_tag = el_node_tags(nd_i)

         if (this%visited(nd_tag) .eq. 0) then  ! once per tag is enough
            this%visited(nd_tag) = el_i

            !nd_i1 = el_node_tags(nd_i)

            tag = p3_tags(nd_i)
            this%v_xy(:, tag) = mesh%v_mshpt_xy(:, nd_tag)
            !this%v_xy(2, tag) = mesh%v_mshpt_xy(2, nd_tag)
            this%v_ety_props(ETY_PROP_PHYSTYPE, tag) = mesh%v_mshpt_physindex(nd_tag)

            !  Vertex => dimension zero
            this%v_ety_props(ETY_PROP_DIMENSION, tag) = 0
         endif
      enddo

      !  now scan the P3 nodes at 1/3 and 2/3 on each edge
      do nd_i=4,6

         nd_tag =el_node_tags(nd_i)

         if(this%visited(nd_tag) .eq. 0) then
            this%visited(nd_tag) = el_i

            ! The P3 edge nodes are 1/3 and 2/3 along the edge
            !  Endpoints of the edge
            k1 = el_node_tags(nd_i-3)
            xx1 = mesh%v_mshpt_xy(1,k1)
            yy1 = mesh%v_mshpt_xy(2,k1)

            k1 = el_node_tags(ip(1,nd_i-3))
            xx2 = mesh%v_mshpt_xy(1,k1)
            yy2 = mesh%v_mshpt_xy(2,k1)

            dx1 = (xx2-xx1) * one_third
            dy1 = (yy2-yy1) * one_third

            !  type of the mid-edge node of the initial P2 mesh

            !  2 nodes per edge (for P3 element)
            do nd_i2=1,2
               k1 = p3_tags(nd_i2+2*(nd_i-4)+3)
               this%v_xy(1,k1) = xx1 + nd_i2*dx1
               this%v_xy(2,k1) = yy1 + nd_i2*dy1
               this%v_ety_props(ETY_PROP_PHYSTYPE,k1) = mesh%v_mshpt_physindex(nd_tag) !  inherit phystype from the mid-edge node of the initial P2 mesh

               this%v_ety_props(ETY_PROP_DIMENSION,k1) = 0      ! these are points, so zero dimension
            enddo
         endif
      enddo


      !  Finally do the interior node

      !  Coordinate of the vertices
      k1 = el_node_tags(1)
      xx1 = mesh%v_mshpt_xy(1,k1)
      yy1 = mesh%v_mshpt_xy(2,k1)

      k1 = el_node_tags(2)
      xx2 = mesh%v_mshpt_xy(1,k1)
      yy2 = mesh%v_mshpt_xy(2,k1)

      k1 = el_node_tags(3)
      xx3 = mesh%v_mshpt_xy(1,k1)
      yy3 = mesh%v_mshpt_xy(2,k1)

      !  Acess the interior node
      k1 = p3_tags(P3_INTERIOR)

      this%v_xy(1, k1) = (xx1+xx2+xx3)*one_third
      this%v_xy(2, k1) = (yy1+yy2+yy3)*one_third

      ! is a point, non-boundary and dimension 0
      this%v_ety_props(ETY_PROP_PHYSTYPE, k1) = 0
      this%v_ety_props(ETY_PROP_DIMENSION, k1) = 0

   enddo

end
