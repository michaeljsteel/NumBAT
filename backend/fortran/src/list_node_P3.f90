


! This Fortran subroutine list_node_P3 is designed to handle the numbering of nodes for P3 (cubic) finite elements in a mesh. Here's a detailed explanation of the subroutine's functionality and structure:

! Purpose:
! The subroutine assigns global numbers to the nodes of P3 elements (which include vertices, edge nodes, and an interior node) for a given mesh. The results are stored in a table that records these node numbers for each element.

! Parameters and Variables:
! Input Parameters:

! nel: Number of elements in the mesh.
! npt: Total number of points (nodes) in the original mesh.
! nnodes: Number of nodes per element.
! n_edge: Number of edges in the mesh.
! Output Parameters:

! npt_p3: Total number of unique nodes after processing P3 elements.
! table_N_E_F: A table that stores node numbers for faces, edges, and the interior node for each element.
! visited: An array used to mark visitedd nodes.
! Internal Variables:

! nnodes_0: Number of vertices in a triangular element.
! nod_el_p: Array to store local node numbers for the current element.
! vert_1, vert_2: Arrays to store vertices of edges for comparison.
! mm, mm2, nn: Counters and offsets for node numbering.
! Description:
! Initialization:

! Initialize the visited array to mark all nodes as unvisitedd.
! Main Loop over Elements:

! Iterate over each element (iel from 1 to nel).
! Vertices Processing:

! For each vertex of the element, check if it has been visitedd.
! If not visitedd, assign a new global number and mark it as visitedd.
! If visitedd, retrieve the global number from the previously processed element.
! Edges Processing:

! For each edge node, determine the vertices of the edge.
! Check if the edge node has been visitedd.
! If not visitedd, assign new global numbers for the edge nodes.
! If visitedd, retrieve the global numbers from the previously processed element, considering the orientation of the edge.
! Interior Node Processing:

! Assign a global number to the single interior node of the P3 element.
! Final Assignment:

! The total number of unique nodes (npt_p3) is updated.
! Structure:
! Initialization of the visited array: This marks all nodes as unvisitedd initially.
! Loop over each element: This processes each element's vertices, edge nodes, and interior node.
! Error Handling: There are checks in place to ensure that nodes and edges are correctly identified and numbered. If an inconsistency is found, the subroutine prints an error message and stops execution.
! Example Use Case:
! This subroutine is typically used in finite element analysis where higher-order elements (like P3 elements) are used to achieve more accurate results. It ensures that each node is uniquely identified across the entire mesh, which is crucial for assembling the global system of equations.

! Notes:
! The subroutine uses integer(8) for large integer values, which allows handling large meshes.
! It uses some Fortran-specific features like parameter and implicit none for better code clarity and safety.
! The subroutine includes detailed error messages and stops execution if a problem is detected, aiding in debugging.
! This subroutine is a critical part of pre-processing in finite element analysis, ensuring that the mesh nodes are correctly numbered and identified for subsequent calculations.


subroutine list_node_P3 (nel, npt, nnodes, n_edge,&
   npt_p3, table_nod, table_N_E_F, visited)


   implicit none
   integer(8) nel, npt, nnodes
   integer(8) n_edge, npt_p3
   integer(8) visited(npt)

   integer(8) table_nod(nnodes,nel), table_N_E_F(14,nel)

!     Local variables
   integer(8) nnodes_0
   parameter (nnodes_0 = 6)
   integer(8) nod_el_p(nnodes_0)
   integer(8) j, k, j1, m1, m2
   integer(8) iel, inod, iel2, inod2, n_face
   integer(8) mm, nn, mm2
   integer(8) vert_1(2), vert_2(2)

   n_face = nel

   visited = 0
   !do j=1,npt  ! OBSELETE
    !  visited(j) = 0
   !enddo

   ! The first 4 entries of table_N_E_F(*,i) correspond to face and edges
   !  These have already been filed in list_face and list_edge

   mm = 4
   mm2 = n_edge + n_face
   nn = 0

   do iel=1,nel

      do inod=1,nnodes
         k = table_nod(inod,iel)
         nod_el_p(inod) = k
      enddo

      !       P3 element: Vertices
      do inod=1,3
         k = nod_el_p(inod)
         if(visited(k) .eq. 0) then
            visited(k) = iel
            nn = nn + 1
            table_N_E_F(inod+mm,iel) = nn + mm2
         else
            iel2 = visited(k)
            inod2 = 0
            do j=1,3
               j1 = table_nod(j,iel2)
               if (k .eq. j1) inod2 = j
            enddo
            if(inod2 .eq. 0) then
               print*, "list_node_P3: problem with a vertex ", iel, inod
               print*, "k, visited(k) = ", k, visited(k)
               stop
            endif
            table_N_E_F(inod+mm,iel) = table_N_E_F(inod2+mm,iel2)
         endif
      enddo

      !       P3 element: nodes on the edges
      do inod=4,6
         k = nod_el_p(inod)
!         Vertices of the edge
         j1 = inod - 3
         vert_1(1) = table_nod(j1,iel)
         j1 = inod - 2
         if (j1 .gt. 3) j1 = j1 - 3
         vert_1(2) = table_nod(j1,iel)
         if(visited(k) .eq. 0) then
            visited(k) = iel
            m1 = 2*(inod-4)+1
            m2 = 2*(inod-4)+2
            do j1=m1,m2
               nn = nn+1
               table_N_E_F(j1+3+mm,iel) = nn + mm2
            enddo
         else
            iel2 = visited(k)
            inod2 = 0

            do j=4,6
               j1=table_nod(j,iel2)
               if (k .eq. j1) then
                  inod2 = j
!               Vertices of the edge
                  j1 = inod2 - 3
                  vert_2(1) = table_nod(j1,iel2)
                  j1 = inod2 - 2
                  if (j1 .gt. 3) j1 = j1 - 3
                  vert_2(2) = table_nod(j1,iel2)
               endif
            enddo

            if(inod2 .eq. 0) then
               print*, "list_node_P3: problem with a node ", iel, inod
               stop
            endif

            do j=1,2
!             local numbering along the edges
               if (vert_2(1) .eq. vert_1(1) .and. &
                  vert_2(2) .eq. vert_1(2)) then
!               The nodes on the edges inod and inod2 are numbered in the same order
!               This is possible only when the elements iel and iel2 have opposite orientations
                  m1 = j+2*(inod-4)+3
                  m2 = j+2*(inod2-4)+3
                  table_N_E_F(m1+mm,iel) = table_N_E_F(m2+mm,iel2)
               elseif (vert_2(1) .eq. vert_1(2) .and. &
                  vert_2(2) .eq. vert_1(1)) then
!               The nodes on the edges inod and inod2 are numbered in the opposite order
!                 ! the numbering of the nodes are reversed
                  j1 = 3 - j
                  m1 = j1+2*(inod-4)+3
                  m2 = j+2*(inod2-4)+3
                  table_N_E_F(m1+mm,iel) = table_N_E_F(m2+mm,iel2)
               else
                  write(*,*) "list_node_P3: problems: ", &
                     "Check the edge endpoints"
                  write(*,*) "inod, table_nod(inod,iel) = ", inod, &
                     table_nod(inod,iel)
                  write(*,*) "inod2, table_nod(inod2,iel2) = ", inod2, &
                     table_nod(inod2,iel2)
                  write(*,*) "iel, iel2 = ", iel, iel2
                  write(*,*) "vert_1 = ", vert_1
                  write(*,*) "vert_2 = ", vert_2
                  write(*,*) "list_node_P3: Aborting..."
                  stop
               endif
            enddo

         endif
      enddo

      !       Numbering the interior nodes of the triangle
!         ! there is only one interior node for a P3 triangle
      do j=1,1   !todo: why is this only one iteration?
         nn = nn+1
         table_N_E_F(j+9+mm,iel) = nn + mm2
      enddo

   enddo

   npt_p3 = nn

   return
end
