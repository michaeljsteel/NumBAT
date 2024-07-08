
! fills out elts table_N_E_F[2..4, 1_n_msh_el]
! each edge is assigned a unique edge number
!
! t_nef[2/3/4,j] = n_msh_el + # of the edge that edge-node is centre of


! table_edge seems to be unused
subroutine label_edges (mesh_props, NEF_props, n_edge, visited)

   use class_MeshProps

   implicit none

   type(MeshProps) :: mesh_props
   type(N_E_F_Props) :: NEF_props

   integer(8) n_edge

   integer(8) visited(mesh_props%n_msh_pts)

   ! ---------------------------------------------

   integer(8) table_edge(4,mesh_props%n_msh_pts)
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

   n_face = mesh_props%n_msh_el
   debug = 0


   visited= 0

   n_edge = 0
   do i=1,mesh_props%n_msh_el

      ! for the edge nodes 4,5,6
      do j=4,mesh_props%nodes_per_el   !  checks some condition on this eleemnt. what is it?
         if (mesh_props%is_boundary_node_2(j,i)) then  ! edge node is on a physical bdy

!            if (type_nod(table_nod(j,i)) .ne. 0) then  ! edge node is on a physical bdy

            ! find node indices (1,2,3) of vertices of the edge
            j1 = list_end(1,j-3)
            j2 = list_end(2,j-3)

            ! Check that both vertices are also bdy points (else would be a broken mesh)
            !if (type_nod(table_nod(j1,i)) .eq. 0 .or. &
            !type_nod(table_nod(j2,i)) .eq. 0) then

            if (.not. mesh_props%is_boundary_node_2(j1,i) .or. &
               .not. mesh_props%is_boundary_node_2(j2,i)) then

               !TODO: hook up error msg
               write(*,*) "list_edge: table_nod = ", &
                  mesh_props%type_node_by_ref(j1,i), &
                  mesh_props%type_node_by_ref(j2,i), &
                  mesh_props%type_node_by_ref(j,i)
               write(*,*) "type_nod(j1) = ", mesh_props%table_nod(j1,i)
               write(*,*) "type_nod(j2) = ", mesh_props%table_nod(j2,i)
               write(*,*) "type_nod(j) = ", mesh_props%table_nod(j,i)
               write(*,*) "list_edge: Aborting..."
               stop
            endif

         endif
      enddo

      !  scan the element edge
      do j=4,mesh_props%nodes_per_el

         j1 = mesh_props%table_nod(j,i)  ! find the node
         k = visited(j1)

         if (k .eq. 0) then        ! a new edge encountered
            n_edge = n_edge + 1

            ! put the three nodes for this edge in table_edge
            visited(j1) = n_edge   ! visited stores its edge number

            j2 = list_end(1,j-3)
            table_edge(1,n_edge) = mesh_props%table_nod(j2,i)

            j2 = list_end(2,j-3)
            table_edge(2,n_edge) = mesh_props%table_nod(j2,i)
            table_edge(3,n_edge) = j1

            !  Table of connectivity for the face (with respect to the triangle element)
            NEF_props%table_nod(j-2,i) = n_edge + n_face
            table_edge(4,n_edge) = n_edge + n_face

         else

            NEF_props%table_nod(j-2,i) = k + n_face
            table_edge(4,k) = k + n_face

         endif

      enddo
   enddo

   return
end
