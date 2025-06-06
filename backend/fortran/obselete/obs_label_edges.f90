
! fills out elts table_N_E_F[2..4, 1_n_msh_elts]
! each edge is assigned a unique edge number
!
! t_nef[2/3/4,j] = n_msh_elts + # of the edge that edge-node is centre of


! table_edge seems to be unused
subroutine label_edges (mesh, m_elnd_to_mshpt, n_edge, visited)

   use class_Mesh

   implicit none

   type(MeshEM) :: mesh

   integer(8) n_edge

   integer(8) visited(mesh%n_msh_pts)
   integer(8) m_elnd_to_mshpt(14,mesh%n_msh_elts)

   ! ---------------------------------------------

   integer(8) table_edge(4,mesh%n_msh_pts)
   integer(8) i, j, k, j1, j2, n_face, debug
   integer(8) list_end(2,3)

   ! ---------------------------------------------

   !  Endpoints of the 3 edges (mid-point) of the reference triangle

   i = 1
   list_end(1,i) = 1
   list_end(2,i) = 2

   i = 2
   list_end(1,i) = 2
   list_end(2,i) = 3

   i = 3
   list_end(1,i) = 1
   list_end(2,i) = 3

   n_face = mesh%n_msh_elts
   debug = 0


   visited= 0

   n_edge = 0
   do i=1,mesh%n_msh_elts

      ! for the edge nodes 4,5,6
      do j=4,mesh%nodes_per_el   !  checks some condition on this eleemnt. what is it?
         if (mesh%is_boundary_mesh_point_by_elt_node(j,i)) then  ! edge node is on a physical bdy

!            if (type_nod(m_elnd_to_mshpt(j,i)) .ne. 0) then  ! edge node is on a physical bdy

            ! find node indices (1,2,3) of vertices of the edge
            j1 = list_end(1,j-3)
            j2 = list_end(2,j-3)

            ! Check that both vertices are also bdy points (else would be a broken mesh)
            !if (type_nod(m_elnd_to_mshpt(j1,i)) .eq. 0 .or. &
            !type_nod(m_elnd_to_mshpt(j2,i)) .eq. 0) then

            if (.not. mesh%is_boundary_mesh_point_by_elt_node(j1,i) .or. .not. mesh%is_boundary_mesh_point_by_elt_node(j2,i)) then

               !TODO: hook up error msg
               write(*,*) "list_edge: m_elnd_to_mshpt = ", &
                  mesh%type_node_by_ref(j1,i), &
                  mesh%type_node_by_ref(j2,i), &
                  mesh%type_node_by_ref(j,i)
               write(*,*) "type_nod(j1) = ", mesh%m_elnd_to_mshpt(j1,i)
               write(*,*) "type_nod(j2) = ", mesh%m_elnd_to_mshpt(j2,i)
               write(*,*) "type_nod(j) = ", mesh%m_elnd_to_mshpt(j,i)
               write(*,*) "list_edge: Aborting..."
               stop
            endif

         endif
      enddo

      !  scan the element edge
      do j=4,mesh%nodes_per_el

         j1 = mesh%m_elnd_to_mshpt(j,i)  ! find the node
         k = visited(j1)

         if (k .eq. 0) then        ! a new edge encountered
            n_edge = n_edge + 1

            ! put the three nodes for this edge in table_edge
            visited(j1) = n_edge   ! visited stores its edge number

            j2 = list_end(1,j-3)
            table_edge(1,n_edge) = mesh%m_elnd_to_mshpt(j2,i)

            j2 = list_end(2,j-3)
            table_edge(2,n_edge) = mesh%m_elnd_to_mshpt(j2,i)
            table_edge(3,n_edge) = j1

            !  Table of connectivity for the face (with respect to the triangle element)
            m_elnd_to_mshpt(j-2,i) = n_edge + n_face
            table_edge(4,n_edge) = n_edge + n_face

         else

            m_elnd_to_mshpt(j-2,i) = k + n_face
            table_edge(4,k) = k + n_face

         endif

      enddo
   enddo

end
