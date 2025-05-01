subroutine MeshEntities_allocate(this, n_msh_el, nberr)

   class(MeshEntities) :: this
   integer(8) :: n_msh_el, n_ddl_ub
   type(NBError) nberr

   ! upper bound for this%n_entities
   n_ddl_ub = 9 * n_msh_el

   this%n_entities = n_ddl_ub  ! provisional value for some memory allocs

   call double_nalloc_2d(this%v_xy, 2_8, n_ddl_ub, 'xy_N_E_F', nberr); RET_ON_NBERR(nberr)
   call integer_nalloc_2d(this%v_ety_props, 2_8, n_ddl_ub, 'type_N_E_F', nberr); RET_ON_NBERR(nberr)
   call integer_nalloc_2d(this%v_tags, 14_8, n_msh_el, 'table_N_E_F', nberr); RET_ON_NBERR(nberr)


   !  Define endpoints of the 3 edges (mid-point) of the reference triangle

   !i = 1
   this%edge_ends(1,1) = 1
   this%edge_ends(2,1) = 2

   !i = 2
   this%edge_ends(1,2) = 2
   this%edge_ends(2,2) = 3

   !i = 3
   this%edge_ends(1,3) = 1
   this%edge_ends(2,3) = 3


end subroutine


 ! !  Storage locations in sequence
 ! !  - tab_N_E_F = table_N_E_F,   shape: 14 x n_msh_el
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

subroutine MeshEntities_build_mesh_tables(this, mesh_raw, nberr)

   use alloc

   class(MeshEntities) :: this
   type(MeshRawEM) :: mesh_raw

   type(NBError) nberr

   ! locals
   integer(8) ui_out
   integer(8), dimension(:), allocatable :: visited

   !  ----------------------------------------------

   ui_out = stdout

   call integer_nalloc_1d(visited, this%n_entities, 'visited', nberr); RET_ON_NBERR(nberr)

   ! Each element has 1 face, 3 edges and 10 P3 nodes
   ! We fill the different rows of v_tags in stages

   !  Fills:  v_tags[1,:]
   call this%count_and_label_faces (mesh_raw%n_msh_el)

   !  Fills: n_edge, v_tags[2:4,:], visited[1:n_msh_pts]
   !         table_edge[1..4,:] (unused)
   call this%count_and_label_edges (mesh_raw, visited, nberr); RET_ON_NBERR(nberr)

   !  Fills: v_tags[5:,:], visited[1:n_msh_pts], n_msh_pts_3
   call this%count_and_label_nodes_P3 (mesh_raw, visited, nberr); RET_ON_NBERR(nberr)

   ! Total number of labelled objects
   this%n_entities = this%n_edges + this%n_faces + this%n_msh_pts_p3

   !  Fills: entities%node_phys_i(1:2, 1:n_ddl), x_E_F(1:2, 1:n_ddl)
   call this%analyse_face_and_edges (mesh_raw, visited)

   call this%analyse_p3_nodes(mesh_raw, visited)

end subroutine


! fills tb_node_labels[1, 1_n_msh_el]
subroutine MeshEntities_count_and_label_faces (this, n_elts)

   class(MeshEntities) :: this
   integer(8) n_elts, i

   ! Every element has one face
   this%n_faces = n_elts

   !  The number of each face is the number of its owning element
   do i=1,n_elts
      this%v_tags(1,i) = i
   enddo

end


! Counts the edges and
! fills tb_node_labels[2..4, 1_n_msh_el]
! corresponding to edges 1,2,3 of each element.
! Each edge is assigned a unique edge number
!    n_msh_el + # of the edge that edge-node is the centre of
! Where an edge lies between two elements it gets the same label
! based on whichever is found first.

! table_edge seems to be unused
subroutine MeshEntities_count_and_label_edges (this, mesh_raw, visited, nberr)

   class(MeshEntities) :: this
   type(MeshRawEM) :: mesh_raw

   integer(8) visited(:)

   type(NBError) nberr

   ! locals
   character(len=EMSG_LENGTH) :: emsg

   integer(8) n_edge
   integer(8) table_edge(4,mesh_raw%n_msh_pts)
   integer(8) iel, jedge
   integer(8) ed_vert_nda, ed_vert_ndb  ! edge vertices nodes
   integer(8) ed_mid_nd ! edge mid node
   integer(8) new_lab, old_lab
   ! ---------------------------------------------

   visited= 0
   n_edge = 0

   ! Check that boundary elements are well constructed.
   ! Not sure why this is an issue
   do iel=1,mesh_raw%n_msh_el  ! for each element

      ! for its edge nodes 4,5,6
      do jedge=4, P2_NODES_PER_EL

         ! if edge node is on a physical bdy
         if (mesh_raw%is_boundary_node_2(jedge,iel)) then

            ! find node indices (1,2,3) of vertices of the edge
            ed_vert_nda = this%edge_ends(1, jedge-3)
            ed_vert_ndb = this%edge_ends(2, jedge-3)

            ! Check that both vertices are also bdy points
            ! (else would be a broken mesh)

            if (.not. mesh_raw%is_boundary_node_2(ed_vert_nda,iel) .or. .not. mesh_raw%is_boundary_node_2(ed_vert_ndb,iel)) then

               write(emsg,*) "list_edge: v_tags = ", &
                  mesh_raw%node_phys_index_by_ref(ed_vert_nda,iel), &
                  mesh_raw%node_phys_index_by_ref(ed_vert_ndb,iel), &
                  mesh_raw%node_phys_index_by_ref(jedge,iel), &
                  "node_phys_i(ed_vert_nda) = ", mesh_raw%elnd_to_mshpt(ed_vert_nda,iel), &
                  "node_phys_i(ed_vert_ndb) = ", mesh_raw%elnd_to_mshpt(ed_vert_ndb,iel), &
                  "node_phys_i(jedge) = ", mesh_raw%elnd_to_mshpt(jedge,iel)
                  call nberr%set(NBERR_BAD_MESH_EDGES, emsg)
               return

            endif

         endif
      enddo

      ! Ok, boundaries are good, let's label the edges

      !  for each edge node 4..6
      do jedge=4,P2_NODES_PER_EL

         ed_mid_nd = mesh_raw%elnd_to_mshpt(jedge,iel)  ! find the node
         old_lab = visited(ed_mid_nd)

         if (old_lab .eq. 0) then        ! a new edge encountered
            n_edge = n_edge + 1

            ! put the three nodes for this edge in table_edge
            new_lab = n_edge+ this%n_faces
            visited(ed_mid_nd) = new_lab  ! visited stores its edge number

            ed_vert_nda = this%edge_ends(1,jedge-3)
            table_edge(1,n_edge) = mesh_raw%elnd_to_mshpt(ed_vert_nda,iel)

            ed_vert_ndb = this%edge_ends(2,jedge-3)
            table_edge(2,n_edge) = mesh_raw%elnd_to_mshpt(ed_vert_ndb,iel)
            table_edge(3,n_edge) = ed_mid_nd
            table_edge(4,n_edge) = new_lab

            !  Table of connectivity for the face (with respect to the triangle element)
            this%v_tags(jedge-2,iel) = new_lab ! -2 because edge 4=edge 1 should end up in row 2

         else  ! Met before, take the
            this%v_tags(jedge-2,iel) = old_lab
            table_edge(4,old_lab- this%n_faces) = old_lab
         endif

      enddo
   enddo

   this%n_edges = n_edge

end

!fills elts
! v_tags[5..7,:] - labels for P3 vertices



subroutine MeshEntities_count_and_label_nodes_P3 (this, mesh_raw, visited, nberr)

   class(MeshEntities) :: this
   type(MeshRawEM) :: mesh_raw

   integer(8) visited(:)

   type(NBError) nberr

   ! Locals
   character(len=EMSG_LENGTH) :: emsg
   integer(8) nod_el_p(P2_NODES_PER_EL)
   integer(8) j, j1, m1, m2
   integer(8) iel, inod, iel2, inod2, nd, nd2
   integer(8) row_off, lab, lab_off
   integer(8) vert_1(2), vert_2(2)

   integer(8) vnd

   visited = 0

   !  The first 4 entries of this%v_tags(*,i) correspond to face and edges
   !  These have already been filed in list_face and list_edge

   row_off = 4                             ! filled rows in v_tags
   lab_off = this%n_edges + this%n_faces  ! assigned labels
   lab = 0

   do iel=1,mesh_raw%n_msh_el

      ! find the absolute node indices of this element
      nod_el_p = mesh_raw%elnd_to_mshpt(:, iel)

      !  P3 element: Vertices
      do inod=1,3
         nd = nod_el_p(inod)          ! find absolute node index
         if(visited(nd) .eq. 0) then  ! if a vew vertex
            visited(nd) = iel         !  tag it with the owning elt

            lab = lab + 1             ! assign a new lab to this vertex
            this%v_tags(inod+row_off,iel) = lab + lab_off
         else
            iel2 = visited(nd)        ! find the owning elt

            ! find which vertex node of the owning elt matches
            inod2 = 0
            do j=1,3
               nd2 = mesh_raw%elnd_to_mshpt(j,iel2)
               if (nd .eq. nd2) inod2 = j
            enddo

            ! illformed mesh
            if(inod2 .eq. 0) then
               write(emsg, *) "list_node_P3: problem with a vertex ", iel, inod, &
                  "nd, visited(nd) = ", nd, visited(nd)
               call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
               return
            endif

            ! copy the vertex label
            this%v_tags(inod+row_off,iel) = this%v_tags(inod2+row_off,iel2)
         endif
      enddo

      !  Have filled 3 more rows, row_off gets +3

      !  P3 element: 6 nodes on the edges surrounding the P2 mid-edge nodes
      do inod=4,6
         nd = nod_el_p(inod)

         !  Find absolute nodes of the vertices of the edge
         vnd = inod - 3
         vert_1(1) = mesh_raw%elnd_to_mshpt(vnd, iel)
         vnd = inod - 2
         if (vnd .gt. 3) vnd = vnd - 3
         vert_1(2) = mesh_raw%elnd_to_mshpt(vnd,iel)

         if (visited(nd) .eq. 0) then  ! new edge
            visited(nd) = iel          ! claim it for this elt
            m1 = 2*(inod-4)+1          ! identify indices 1.2, 3.4, or 5.6
            m2 = 2*(inod-4)+2
            do j1=m1,m2
               lab = lab+1
               this%v_tags(j1+3+row_off,iel) = lab + lab_off
            enddo
         else
            iel2 = visited(nd)        ! find the owning elt
            inod2 = 0

            ! find which edge node of the owning elt matches
            do j=4,6
               nd2=mesh_raw%elnd_to_mshpt(j,iel2)

               if (nd .eq. nd2) then
                  inod2 = j
                  !  Vertices of the edge
                  vnd = inod2 - 3
                  vert_2(1) = mesh_raw%elnd_to_mshpt(vnd, iel2)
                  vnd = inod2 - 2
                  if (vnd .gt. 3) vnd = vnd - 3
                  vert_2(2) = mesh_raw%elnd_to_mshpt(vnd, iel2)
               endif

            enddo

            if(inod2 .eq. 0) then
               write(emsg, *)  "list_node_P3: problem with a node ", iel, inod
               call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
            endif

            ! copy the labels
            do j=1,2
               !  local numbering along the edges
               if (vert_2(1) .eq. vert_1(1) .and. vert_2(2) .eq. vert_1(2)) then
                  !  The nodes on the edges inod and inod2 are numbered in the same order
                  !  This is possible only when the elements iel and iel2 have opposite orientations
                  m1 = j+2*(inod-4)
                  m2 = j+2*(inod2-4)
               elseif (vert_2(1) .eq. vert_1(2) .and.  vert_2(2) .eq. vert_1(1)) then
                  !  The nodes on the edges inod and inod2 are numbered in the opposite order
                  !  the numbering of the nodes are reversed
                  j1 = 3 - j
                  m1 = j1+2*(inod-4)
                  m2 = j+2*(inod2-4)
               else
                  write(emsg,*) "list_node_P3: problems: ", &
                     "Check the edge endpoints", &
                     "inod, elnd_to_mshpt(inod,iel) = ", inod, &
                     mesh_raw%elnd_to_mshpt(inod,iel), &
                     "inod2, elnd_to_mshpt(inod2,iel2) = ", inod2, &
                     mesh_raw%elnd_to_mshpt(inod2,iel2), &
                     "iel, iel2 = ", iel, iel2, &
                     "vert_1 = ", vert_1, &
                     "vert_2 = ", vert_2
                  call nberr%set(NBERR_BAD_MESH_VERTICES, emsg)
                  return
               endif

               this%v_tags(m1+row_off+3,iel) = this%v_tags(m2+row_off+3,iel2)
            enddo

         endif
      enddo

      !  Have filled 6 more rows, row_off gets +9

      !  Numbering the interior nodes of the triangle
      !  there is only one interior node for a P3 triangle
      lab = lab+1
      this%v_tags(1+9+row_off,iel) = lab + lab_off

   enddo

   this%n_msh_pts_p3 = lab

end



subroutine MeshEntities_analyse_face_and_edges (this, mesh_raw, visited)


   class(MeshEntities) :: this
   type(MeshRawEM) :: mesh_raw

   integer(8) visited(:)

   double precision, parameter :: one_third = 1.d0/3.d0
   integer(8) iel, j, nd, tag
   integer(8) type_n(10)
   double precision el_xy(2,6)

   this%v_ety_props = 0
   visited= 0


   do iel=1,mesh_raw%n_msh_el

      ! find the elt's material and its node's locations
      do j=1,P2_NODES_PER_EL
         nd = mesh_raw%elnd_to_mshpt(j,iel)
         type_n(j) = mesh_raw%v_nd_physindex(nd)
         el_xy(:,j) = this%v_xy(:,nd)
      enddo

      ! Face properties
      tag = this%v_tags(1, iel)
      ! Position is the barycentre
      this%v_xy(:,tag) = (el_xy(:,1) + el_xy(:,2) + el_xy(:,3)) * one_third

      this%v_ety_props(1,tag) = 0  !  Topologically, a face is an interior domain
      this%v_ety_props(2,tag) = 2  !  Face => dimension two


      !  scan the 3 element edges
      do j=1,3
         tag = this%v_tags(j+1,iel)  ! edges start at row 2

         this%v_xy(:,tag) = el_xy(:,j+3)

         if (visited(tag) .eq. 0) then  ! only do each tag once
            visited(tag) = 1
            this%v_ety_props(1,tag) = type_n(j+3)

            !  Edge => dimension one
            this%v_ety_props(2,tag) = 1
         endif

      enddo
   enddo
end




subroutine MeshEntities_analyse_p3_nodes(this, mesh_raw, visited)

   class(MeshEntities) :: this
   type(MeshRawEM) :: mesh_raw

   integer(8) visited(:)


   integer(8)  k1, n, ind, ip(2,3), nd, tag
   integer(8) iel, inod,  inod2, row_off
   integer(8) el_nodes(6), p3_tags(N_ENTITY_PER_EL)

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


   visited = 0


   !  The first 4 entries of this%v_tags(*,i) correspond to face and P2 edges and have been done
   row_off = 4
   do iel=1,mesh_raw%n_msh_el

      ! do inod=1,P2_NODES_PER_EL
      !    el_nodes(inod) = mesh_raw%elnd_to_mshpt(inod,iel)
      ! enddo
      el_nodes = mesh_raw%elnd_to_mshpt(:, iel)

      !  the 10 node of a P3 element
      ! do inod=1,P3_NODES_PER_EL
      !    p3_tags(inod) = this%v_tags(inod+row_off,iel)
      ! enddo
      p3_tags(1:P3_NODES_PER_EL)  = this%v_tags(row_off+1:row_off+10,iel)


      !  scan the P3 vertices
      do inod=1,3
         nd = el_nodes(inod)

         if (visited(nd) .eq. 0) then  ! once per tag is enough
            visited(nd) = iel
            !inod1 = el_nodes(inod)
            tag = p3_tags(inod)
            this%v_xy(1, tag) = mesh_raw%v_nd_xy(1, nd)
            this%v_xy(2, tag) = mesh_raw%v_nd_xy(2, nd)
            this%v_ety_props(1, tag) = mesh_raw%v_nd_physindex(nd)

            !  Vertex => dimension zero
            this%v_ety_props(2, tag) = 0
         endif
      enddo

      !  scan the P3 nodes located on the edges, tags 4..9 + row_off
      do inod=4,6

         nd =el_nodes(inod)

         if(visited(nd) .eq. 0) then
            visited(nd) = iel

            ! The P3 edge nodes are 1/3 and 2/3 along the edge
            !  Endpoints of the edge
            k1 = el_nodes(inod-3)
            xx1 = mesh_raw%v_nd_xy(1,k1)
            yy1 = mesh_raw%v_nd_xy(2,k1)

            k1 = el_nodes(ip(1,inod-3))
            xx2 = mesh_raw%v_nd_xy(1,k1)
            yy2 = mesh_raw%v_nd_xy(2,k1)

            dx1 = (xx2-xx1) * one_third
            dy1 = (yy2-yy1) * one_third

            !  type of the mid-edge node of the initial P2 mesh
            ind = mesh_raw%v_nd_physindex(nd)

            !  2 nodes per edge (for P3 element)
            do inod2=1,2
               k1 = p3_tags(inod2+2*(inod-4)+3)
               this%v_xy(1,k1) = xx1 + inod2*dx1
               this%v_xy(2,k1) = yy1 + inod2*dy1
               this%v_ety_props(1,k1) = ind

               !  Node => dimension zero
               this%v_ety_props(2,k1) = 0
            enddo
         endif
      enddo

      !  Coordinate of the vertices
      k1 = el_nodes(1)
      xx1 = mesh_raw%v_nd_xy(1,k1)
      yy1 = mesh_raw%v_nd_xy(2,k1)

      k1 = el_nodes(2)
      xx2 = mesh_raw%v_nd_xy(1,k1)
      yy2 = mesh_raw%v_nd_xy(2,k1)

      k1 = el_nodes(3)
      xx3 = mesh_raw%v_nd_xy(1,k1)
      yy3 = mesh_raw%v_nd_xy(2,k1)

      !  The tenth node is at the center of the triangle
      !  dimension(P3) = 10
      n = 10

      k1 = p3_tags(n)

      this%v_xy(1,k1) = (xx1+xx2+xx3)*one_third
      this%v_xy(2,k1) = (yy1+yy2+yy3)*one_third

      !  interior node
      this%v_ety_props(1,k1) = 0

      !  Node => dimension zero
      this%v_ety_props(2,k1) = 0

   enddo

end




subroutine MeshEntitiesAC_allocate(this, n_msh_el, nberr)

   class(MeshEntitiesAC) :: this
   integer(8) :: n_msh_el, n_ddl_ub
   type(NBError) nberr

   ! upper bound for this%n_entities
   n_ddl_ub = 9 * n_msh_el

   this%n_entities = n_ddl_ub  ! provisional value for some memory allocs

   call double_nalloc_2d(this%v_xy, 2_8, n_ddl_ub, 'xy_N_E_F', nberr); RET_ON_NBERR(nberr)
   call integer_nalloc_2d(this%v_ety_props, 2_8, n_ddl_ub, 'type_N_E_F', nberr); RET_ON_NBERR(nberr)
   call integer_nalloc_2d(this%v_tags, 14_8, n_msh_el, 'table_N_E_F', nberr); RET_ON_NBERR(nberr)


   !  Define endpoints of the 3 edges (mid-point) of the reference triangle

   !i = 1
   this%edge_ends(1,1) = 1
   this%edge_ends(2,1) = 2

   !i = 2
   this%edge_ends(1,2) = 2
   this%edge_ends(2,2) = 3

   !i = 3
   this%edge_ends(1,3) = 1
   this%edge_ends(2,3) = 3


end subroutine
