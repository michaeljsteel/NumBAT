#include "numbat_decl.h"

!  This file performs RCM reordering to make the bandwidth of the adjacency matrix smaller.
!  This makes sparse matrix operations more efficient.
!  The shuffling is based on the triangular elements in v_gmsh_tri_nodes only.
!  But we must also permute the indices and positions in v_nd_iphyscurve and vx, vy.

subroutine balance_fem_node_graph(n_pts, n_gmsh_tri, v_gmsh_tri_nodes, v_nd_iphyscurve, &
   vx, vy, assertions_on, errco, emsg)

   use numbatmod

   integer assertions_on, errco
   character emsg*(STRINGLEN)
   
   integer n_pts, n_gmsh_tri
   integer v_gmsh_tri_nodes(6,n_gmsh_tri), v_nd_iphyscurve(n_pts)
   double precision vx(n_pts), vy(n_pts)

   integer num_adj               ! total accumulated number of adj nodes
   integer xadj(MAX_N_PTS+1)     ! cumulative number of connections
   integer adjncy(MAX_LONG_ADJ)  ! how many nodes have a connection to this one

   ! builds xadj and num_adj
   call make_adjacency_cumulative_vector(n_pts, n_gmsh_tri, v_gmsh_tri_nodes, xadj, num_adj)

   if (assertions_on .ne. 0) then
   call assert_no_larger_than(num_adj, MAX_LONG_ADJ, 'renumber_nodes','num_adj <= MAXLONGADJ', -11, errco, emsg)
   RETONERROR(errco)
   endif

   ! builds adjncy
   call make_adjacency_matrix(n_pts,n_gmsh_tri, v_gmsh_tri_nodes, xadj, num_adj, adjncy)
!         write(*,*) 'renumbering with: num_adj = ', num_adj

   call rebalance_adjacency_matrix(n_pts, n_gmsh_tri, v_gmsh_tri_nodes, xadj, num_adj, &
      adjncy, vx, vy, v_nd_iphyscurve, errco, emsg)

   return
end

!-------------------------------------------------------------------------------
! builds xadj and num_adj from v_gmsh_tri_nodes
subroutine make_adjacency_cumulative_vector(n_pts, n_gmsh_tri, v_gmsh_tri_nodes, xadj, num_adj)

   use numbatmod

   integer n_gmsh_tri, n_pts
   integer v_gmsh_tri_nodes(6,n_gmsh_tri)
   integer xadj(MAX_N_PTS+1)     ! cumulative number of connections
   integer num_adj

   integer i, j, nd, i1
   integer t_nodes(6)

   integer visited(MAX_N_PTS), neighbours(MAX_N_PTS)

   do i = 1,n_pts
      visited(i) = 0
      neighbours(i)=0
   enddo

   ! Find the number of associations each triangle node has

   do j = 1,n_gmsh_tri     ! for all triangle (6 node) elements

      do i=1,6              ! copy the 6 nodes for this element
         t_nodes(i) = v_gmsh_tri_nodes(i,j)
         nd = t_nodes(i)
         neighbours(nd) = neighbours(nd) + 6 - 1    ! account for the 5 other nodes of this element
      enddo

      do i=1,3
         i1 = t_nodes(i+3)            ! consider the edge nodes
         if(visited(i1) .eq. 0) then  !first time visiting this edge node
            visited(i1) = 1
         else
            ! So the points on edge i have already been counted, and must be un-double counted
            ! The three points on the edge have already met their two neighbours on this edge, so we subtract 2.
            ! Nodes in order clockwise around the triangle are 1(v) 4(e) 2(v) 5(e) 3(v) 6(e)
            nd = t_nodes(i)         ! the first vertex
            neighbours(nd) = neighbours(nd) - 2
            nd = t_nodes(1+mod(i,3))          ! the second vertex, wrapping back to 1 if it is vertex 3
            neighbours(nd) = neighbours(nd) - 2
            nd = t_nodes(i+3)       ! the edge
            neighbours(nd) = neighbours(nd) - 2
         endif
      enddo
   enddo

   !do i = 1,n_pts
   !   write(*,*) 'neighbours', i, neighbours(i)
   !enddo

   ! make xadj the vector of cumulative sum of neighbours
   xadj(1) = 1
   do i = 1,n_pts
      xadj(i+1) = xadj(i) + neighbours(i)
   enddo
   num_adj = xadj(n_pts+1)

   end subroutine

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Make adjncy from xadj, num_adj, v_gmsh_tri_nodes

subroutine make_adjacency_matrix(n_pts, n_elts_tri,  v_gmsh_tri_nodes, xadj, num_adj, adjncy)

   use numbatmod
   integer n_elts_tri, n_pts, num_adj
   integer v_gmsh_tri_nodes(6,n_elts_tri)
   integer xadj(n_pts+1), adjncy(num_adj)

   integer lb2(MAX_N_PTS)
   integer t_nodes(6)
   integer i, j, k, i1, k1, ind1, ind2, m, m1

   do i = 1,n_pts
      lb2(i)=0
   enddo

   do j = 1,n_elts_tri     ! For all (6 node) triangles
      do i=1,6
         t_nodes(i) = v_gmsh_tri_nodes(i,j)  ! get the current 6 nodes of this elt
      enddo

      do i=1,6              ! for each of the 6 nodes
         i1=t_nodes(i)      !  find its node index in 1..n_pts
         ind1=xadj(i1)      !  find its cumulative adjacencies

         do k=1,6           !    comparing it against all 6 nodes at this elt
            k1=t_nodes(k)
            ind2=0               !  if self compare, ind2=1
            if(k .eq. i) then
               ind2=1
               goto 5
            endif

            do m=1,lb2(i1)
               m1=adjncy(ind1+m-1)
               if(m1 .eq. k1) then
                  ind2=1
                  goto 5
               endif
            enddo

5           continue

            if(ind2 .eq. 0) then
               lb2(i1) = lb2(i1)+1
               adjncy(ind1+lb2(i1)-1)=k1
            endif
         enddo
      enddo
   enddo

end subroutine

!----------------------------------------------------------------------------
! Every elt of adjncy except the last should be nonzero (apparently).
subroutine assert_good_adjmat_permutation(num_adj, adjncy, errco, emsg)
   use numbatmod
   integer errco
   character emsg*(STRINGLEN)
   integer num_adj, i
   integer adjncy(num_adj)

   do i = 1,num_adj-1
      if(adjncy(i) .eq. 0) then
         write(emsg,*) 'rebalance_adjaceny_matrix ',  '(in conv_gmsh.f): ', &
            ' Attention, permuted adjacency matrix has unexpected zero values', 'i, adjncy(i) = ', i, adjncy(i)
         errco=-53
         return
      endif
   enddo
   end subroutine


!----------------------------------------------------------------------------
! Shuffle v_gmsh_tri_nodes, vx, vy 
subroutine apply_node_permutations(n_pts, n_elts_tri, perm, idfn, v_gmsh_tri_nodes, vx, vy)
   use numbatmod

   integer n_elts_tri, n_pts
   integer perm(MAX_N_PTS)
   integer idfn(n_pts)
   integer v_gmsh_tri_nodes(6,n_elts_tri)
   double precision vx(n_pts), vy(n_pts)

   integer i,j
   integer invperm(MAX_N_PTS), idfn_r(MAX_N_PTS)
   double precision x_r(MAX_N_PTS), y_r(MAX_N_PTS)
   
      do i = 1, n_pts
         invperm(perm(i)) = i
      enddo
   
      do i = 1, n_elts_tri
         do j = 1, 6
            v_gmsh_tri_nodes(j,i) = invperm(v_gmsh_tri_nodes(j,i))
         enddo
      enddo
   
      do i = 1, n_pts
         idfn_r(i) = idfn(perm(i))
         x_r(i) = vx(perm(i))
         y_r(i) = vy(perm(i))
      enddo
   
      do i = 1, n_pts
         idfn(i) = idfn_r(i)
         vx(i) = x_r(i)
         vy(i) = y_r(i)
      enddo
   
   end subroutine


!c###############################################
! Builds adjacency vector adjcny from xadj
! Shuffles vx, vy, idfn accordingly

! This Fortran subroutine, named `rebalance_adjacency_matrix`, appears to perform the following tasks:

! 1. **Building the adjacency matrix for a mesh of triangular elements:**
!    - It takes as input the number of points (`n_pts`), the number of triangular elements (`n_elts_tri`), the connectivity array of nodes for each element (`v_gmsh_tri_nodes`), the node identification array (`idfn`), the cumulative sum of degrees array (`xadj`), and the adjacency list (`adjncy`).
!    - For each triangular element, it updates the adjacency list (`adjncy`) based on the connectivity of nodes. It avoids duplicates and self-comparisons.

! 2. **Applying the Reverse Cuthill-McKee (RCM) algorithm to renumber nodes:**
!    - It initializes arrays (`mask`, `perm`, `xls`, etc.) for use in the RCM algorithm.
!    - It calls the `genrcm` subroutine, which is responsible for applying the RCM algorithm to generate a new ordering of nodes (`perm`).
!    - It checks for potential issues in the adjacency list (`adjncy`) and returns an error if necessary.

! 3. **Reordering mesh-related arrays based on the RCM permutation:**
!    - It creates an inverse permutation array (`invperm`) to map the new ordering back to the original ordering.
!    - It applies the inverse permutation to the connectivity array of nodes for each triangular element (`v_gmsh_tri_nodes`).
!    - It applies the inverse permutation to the node identification array, `vx` (x-coordinates), and `vy` (y-coordinates).

! 4. **Additional notes:**
!    - The subroutine includes error checking for potential issues with the adjacency list.
!    - It provides the option to reorganize additional arrays based on the RCM ordering (such as `idfn`, `vx`, and `vy`), but this is conditional on the logical value `.true.`.

! It's worth noting that the RCM algorithm is used here to reduce the bandwidth of the adjacency matrix, which can be beneficial for certain numerical computations involving sparse matrices. The reordering of arrays based on the RCM permutation is a common step in preparing data for efficient sparse matrix computations.


subroutine rebalance_adjacency_matrix(n_pts, n_elts_tri, v_gmsh_tri_nodes, xadj, num_adj, adjncy, vx, vy, idfn,  errco, emsg)

   use numbatmod

   integer errco
   character emsg*(STRINGLEN)

   integer n_elts_tri, n_pts, num_adj
   integer v_gmsh_tri_nodes(6,n_elts_tri), idfn(n_pts)
   integer xadj(n_pts+1), adjncy(num_adj)
   double precision vx(n_pts), vy(n_pts)
   

   integer mask(MAX_N_PTS), perm(MAX_N_PTS)
   integer xls(MAX_N_PTS)
   

   call genrcm(n_pts, xadj, adjncy, perm, mask, xls)


   call assert_good_adjmat_permutation(num_adj, adjncy, errco, emsg)
   RETONERROR(errco)
   
   call apply_node_permutations(n_pts, n_elts_tri, perm, idfn, v_gmsh_tri_nodes, vx, vy)

end subroutine



!--------------------------------------------------------------------------------------------------------

! This Fortran subroutine, named `genrcm`, seems to be generating a Reverse Cuthill-McKee (RCM) ordering 
! for a graph represented by its adjacency list. The RCM algorithm is used to reduce the bandwidth of a 
! sparse matrix, which can be beneficial for solving sparse linear systems more efficiently.

! Here's a breakdown of the key components of the function:

! - `neqns`: Input parameter representing the number of equations (nodes) in the graph.
! - `xadj`: Cumulative sum of the degrees of nodes in the graph.
! - `adjncy`: Adjacency list representation of the graph.
! - `perm`: Output parameter representing the RCM ordering of the nodes.
! - `mask`: An array indicating whether a node has been visited (`mask(i) = 1` means not visited, and `mask(i) = 0` means visited).
! - `xls`: An array storing the cumulative sum of nodes at each level in the level structure.

! The subroutine initializes a mask array where all nodes are marked as not visited (`mask(i) = 1`). 
! It then iterates over the nodes (`neqns`) to find the RCM ordering. For each unvisited node, it 
! invokes the `fnroot` subroutine to generate a level structure (`xls`) using the breadth-first search (BFS) 
! approach and updates the mask accordingly. 
! Subsequently, the `rcm` subroutine is called to perform the Reverse Cuthill-McKee ordering starting from the 
! specified root node.

! The ordering obtained from the `rcm` subroutine is stored in the `perm` array. The process continues until all nodes are visited. 
! The loop is designed to handle disconnected components in the graph.

! In summary, this subroutine applies the Reverse Cuthill-McKee algorithm to generate a node ordering (`perm`) 
! that minimizes the bandwidth of the graph represented by the adjacency list. 
! The `xls` array is used to store the level structure, and the `mask` array keeps track of visited nodes during the process.

subroutine  genrcm (neqns, xadj, adjncy, perm, mask, xls)

      integer adjncy(*), mask(*), perm(*), xls(*)
      integer xadj(*), ccsize, i, neqns, nlvl, num, root
   
      do i = 1, neqns
         mask(i) = 1
      enddo
      
      num = 1
      do i = 1, neqns
         if (mask(i).eq.0) continue

         root = i
         call fnroot (root,xadj,adjncy,mask,nlvl,xls,perm(num))
         call rcm (root,xadj,adjncy,mask,perm(num),ccsize,xls)

         num = num + ccsize
         if (num.gt.neqns) return

      enddo
   end
   
!    ------------------------------------------------------------------

subroutine fnroot (root,xadj,adjncy,mask,nlvl,xls,ls)

! This Fortran subroutine appears to be part of a larger algorithm for finding a root node for a level structure in a graph. The subroutine aims to identify an appropriate root node that minimizes the degree of its neighbors, possibly optimizing for a more balanced level structure. Let's break down the key components:

! - `root`: The input and output parameter representing the current root node.
! - `xadj`: Cumulative sum of the degrees of nodes in the graph.
! - `adjncy`: Adjacency list representation of the graph.
! - `mask`: An array indicating whether a node has been visited (`mask(i) > 0` means visited, and `mask(i) = 0` means not visited).
! - `nlvl`: Input parameter representing the current level of the structure.
! - `xls`: An array storing the cumulative sum of the nodes at each level in the level structure.
! - `ls`: An array representing the level structure.

! The subroutine begins by calling another subroutine `rootls` (not provided here), which appears to be responsible for generating the level structure starting from the specified root node. The resulting level structure is stored in `xls` and `ls`.

! After generating the initial level structure, the subroutine evaluates whether the current level (`nlvl`) is either 1 or equal to the size of the connected component (`ccsize`). If so, it returns without further processing, assuming that the level structure is already optimal.

! If the current level is not the first or equal to the connected component size, the subroutine enters a loop (`100`) where it iterates over the nodes in the current level and evaluates their degrees and the degrees of their neighbors. It selects a new root node based on the node with the minimum degree among its neighbors.

! The process continues until it iterates through all nodes in the current level. After identifying a new root node, it calls `rootls` again with the updated root node to generate a new level structure (`nunlvl`). If the new level structure has fewer levels than the current level (`nunlvl <= nlvl`), the subroutine returns. Otherwise, it updates the current level (`nlvl`) and repeats the process until the level structure is considered optimal.

! In summary, this subroutine aims to find an optimal root node for a level structure in a graph, considering the degrees of nodes and their neighbors. The goal seems to be achieving a more balanced level structure in the graph.



!    ------------------------------------------------------------------

   integer adjncy(*), ls(*), mask(*), xls(*)
   integer xadj(*), ccsize, j, jstrt, k, kstop, kstrt
   integer mindeg, nabor, ndeg, nlvl, node, nunlvl, root

!    ------------------------------------------------------------------

   call rootls (root,xadj,adjncy,mask,nlvl,xls,ls)
   ccsize = xls(nlvl+1) - 1
   if (nlvl.eq.1 .or. nlvl.eq.ccsize) return
100 jstrt = xls(nlvl)
   mindeg = ccsize
   root = ls(jstrt)
   if (ccsize.eq.jstrt) go to 400
   do 300 j = jstrt, ccsize
      node = ls(j)
      ndeg = 0
      kstrt = xadj(node)
      kstop = xadj(node+1) - 1
      do 200 k = kstrt, kstop
         nabor = adjncy(k)
         if (mask(nabor).gt.0) ndeg = ndeg + 1
200   continue
      if (ndeg.ge.mindeg) go to 300
      root = node
      mindeg = ndeg
300 continue
400 call rootls (root,xadj,adjncy,mask,nunlvl,xls,ls)
   if (nunlvl.le.nlvl) return
   nlvl = nunlvl
   if (nlvl.lt.ccsize) go to 100


end

!c###############################################


subroutine  rcm (root,xadj,adjncy,mask,perm,ccsize,deg)

!         This Fortran subroutine, named `rcm`, implements the Reverse Cuthill-McKee (RCM) algorithm. The RCM algorithm is used for reordering the nodes of a graph to reduce the bandwidth of its adjacency matrix. Lower bandwidth can be beneficial for certain numerical computations, particularly in the context of sparse matrix operations.

! Here's a breakdown of the key components of the `rcm` subroutine:

! - `root`: Input parameter representing the starting node for the RCM traversal.
! - `xadj`: Cumulative sum of the degrees of nodes in the graph.
! - `adjncy`: Adjacency list representation of the graph.
! - `mask`: An array indicating whether a node has been visited (`mask(i) = 0` means visited, and `mask(i) = 1` means not visited).
! - `perm`: Output parameter representing the RCM ordering of the nodes.
! - `ccsize`: Output parameter representing the size of the connected component.
! - `deg`: Output parameter representing the degrees of nodes in the connected component.

! The subroutine proceeds as follows:

! 1. It calls another subroutine named `degree` to calculate the degrees of nodes in the connected component starting from the specified root. The results are stored in the `deg` array, and the connected component size is returned in `ccsize`.

! 2. It initializes the mask for the root node and checks if the connected component size is less than or equal to 1. If so, there is no reordering needed, and the subroutine returns.

! 3. It sets up a loop (`100`) to perform the RCM traversal. The traversal is essentially a breadth-first search (BFS) on the graph.

! 4. Within the loop:
!    - It iterates over the nodes in the current level and updates the mask, perm, and lnbr (last node in the current level).
!    - It then reorders the nodes in the current level based on their degrees using a modified insertion sort.

! 5. The loop continues until all nodes in the connected component are processed.

! 6. Finally, it performs a final reordering to ensure symmetry in the permutation.

! In summary, this subroutine applies the Reverse Cuthill-McKee algorithm to generate a reordering (`perm`) that reduces the bandwidth of the graph's adjacency matrix. The degree information is used to influence the ordering during the BFS traversal. The output parameters include the size of the connected component (`ccsize`) and the degrees of nodes in the connected component (`deg`).
!    ------------------------------------------------------------------

   integer adjncy(*),deg(*),mask(*),perm(*)
   integer xadj(*),ccsize,fnbr,i,j,jstop,jstrt,k,l,lbegin,lnbr
   integer lperm,lvlend,nbr,node,root

!    ------------------------------------------------------------------

   call degree (root,xadj,adjncy,mask,deg,ccsize,perm)
   mask(root) = 0
   if (ccsize.le.1) return
   lvlend = 0
   lnbr = 1
100 lbegin = lvlend + 1
   lvlend = lnbr
   do 600 i = lbegin, lvlend
      node = perm(i)
      jstrt = xadj(node)
      jstop = xadj(node+1) - 1
      fnbr = lnbr + 1
      do 200 j = jstrt, jstop
         nbr = adjncy(j)
         if (mask(nbr).eq.0) go to 200
         lnbr = lnbr + 1
         mask(nbr) = 0
         perm(lnbr) = nbr
200   continue
      if (fnbr.ge.lnbr) go to 600
      k = fnbr
300   l = k
      k = k + 1
      nbr = perm(k)
400   if (l.lt.fnbr) go to 500
      lperm = perm(l)
      if (deg(lperm).le.deg(nbr)) go to 500
      perm(l+1) = lperm
      l = l - 1
      go to 400
500   perm(l+1) = nbr
      if (k.lt.lnbr) go to 300
600 continue
   if (lnbr.gt.lvlend) go to 100
   k = ccsize/2
   l = ccsize
   do 700 i = 1, k
      lperm = perm(l)
      perm(l) = perm(i)
      perm(i) = lperm
      l = l - 1
700 continue


end
!c###############################################
!    ------------------------------------------------------------------

subroutine  rootls (root,xadj,adjncy,mask,nlvl,xls,ls)

!         This Fortran subroutine `rootls` appears to perform a level structure traversal of a graph starting from a specified root node. It uses a breadth-first search (BFS) approach to explore the graph and organize nodes into levels.

! Here's a breakdown of the key components of the function:

! - `root`: Input parameter representing the starting node for the BFS traversal.
! - `xadj`: Cumulative sum of the degrees of nodes in the graph.
! - `adjncy`: Adjacency list representation of the graph.
! - `mask`: An array indicating whether a node has been visited (`mask(i) = 0` means not visited, and `mask(i) = 1` means visited).
! - `nlvl`: Output parameter representing the total number of levels in the level structure.
! - `xls`: An array storing the cumulative sum of nodes at each level in the level structure.
! - `ls`: An array representing the level structure.

! The subroutine initializes by marking the root node as visited (setting `mask(root) = 0`), and adding it to the level structure (`ls` array). It then enters a loop (`200`) that iteratively explores the neighbors of the nodes in the current level. The loop continues until all nodes at the current level have been processed.

! Within the loop:
! 1. The starting and ending indices for the current level are determined (`lbegin` and `lvlend`).
! 2. The neighbors of each node in the current level are explored, and if a neighbor has not been visited (`mask(nbr) = 1`), it is added to the level structure.
! 3. The `mask` array is updated to mark visited nodes.

! After processing all nodes in the current level, the subroutine checks if there are more nodes to explore in subsequent levels. If so, it increments the level count (`nlvl`), updates the cumulative sum array (`xls`), and repeats the process.

! The subroutine continues until all nodes in the graph are visited. The final result is a level structure stored in the `ls` array, with information about the cumulative sum of nodes at each level stored in the `xls` array. The `mask` array is used to keep track of visited nodes during the traversal.

!    ------------------------------------------------------------------

   integer adjncy(*), ls(*), mask(*), xls(*)
   integer xadj(*), i, j, jstop, jstrt, lbegin
   integer ccsize, lvlend, lvsize, nbr, nlvl, node, root
!    ------------------------------------------------------------------

   mask(root) = 0
   ls(1) = root
   nlvl = 0
   lvlend = 0
   ccsize = 1
200 lbegin = lvlend + 1
   lvlend = ccsize
   nlvl = nlvl + 1
   xls(nlvl) = lbegin
   do 400 i = lbegin, lvlend
      node = ls(i)
      jstrt = xadj(node)
      jstop = xadj(node + 1) - 1
      if (jstop.lt.jstrt) go to 400
      do 300 j = jstrt, jstop
         nbr = adjncy(j)
         if (mask(nbr).eq. 0) go to 300
         ccsize = ccsize + 1
         ls(ccsize) = nbr
         mask(nbr) = 0
300   continue
400 continue
   lvsize = ccsize - lvlend
   if (lvsize.gt.0) go to 200
   xls(nlvl+1) = lvlend + 1
   do 500 i = 1, ccsize
      node = ls(i)
      mask(node) = 1
500 continue


end


subroutine  degree (root,xadj,adjncy,mask,deg,ccsize,ls)

   !         ChatGPT says:
   !         This Fortran subroutine appears to be implementing a breadth-first search (BFS) traversal of a graph starting from a given `ROOT` node. The graph is represented using an adjacency list with arrays `XADJ` and `ADJNCY`, where `XADJ` contains the cumulative sum of the degrees of nodes, and `ADJNCY` contains the adjacency information for each node.
   
   ! Here's a breakdown of the key components of the function:
   
   ! - `ROOT`: The starting node for the BFS traversal.
   ! - `XADJ`: Cumulative sum of the degrees of nodes in the graph.
   ! - `ADJNCY`: Adjacency list representation of the graph.
   ! - `MASK`: An array indicating whether a node has been visited (`MASK(i)=0` means not visited, and `MASK(i)=1` means visited).
   ! - `DEG`: An array to store the degree of each node.
   ! - `CCSIZE`: A variable to keep track of the size of the connected component.
   ! - `LS`: An array representing the BFS traversal order.
   
   ! The function initializes with the `ROOT` node and performs BFS by traversing the graph layer by layer. The BFS traversal is done in a loop (`100`) until all nodes are visited. The process involves updating the degree of each node (`DEG` array), marking visited nodes in the `MASK` array, and updating the BFS traversal order in the `LS` array.
   
   ! The key steps within the loop (`100`) are:
   
   ! 1. Setting the beginning and end indices for the current layer of the BFS traversal.
   ! 2. Iterating over the nodes in the current layer.
   ! 3. For each node, updating its degree (`IDEG`), marking it as visited, and updating the BFS traversal order.
   
   ! The loop continues until all nodes in the current layer are processed. The process then repeats for the next layer until the entire graph is traversed.
   
   ! Finally, the function resets the signs of the `XADJ` array for all nodes, effectively undoing the changes made during the traversal.
   
   ! It's worth noting that this subroutine does not explicitly handle disconnected components in the graph. It assumes that the graph is connected. If there are multiple connected components, you might need to invoke this subroutine for each component separately.
   
   
      integer adjncy(*), deg(*), ls(*), mask(*)
      integer xadj(*), ccsize, i, ideg, j, jstop, jstrt
      integer lbegin, lvlend, lvsize, nbr, node, root
   
   
      ls(1) = root
      xadj(root) = -xadj(root)
      lvlend = 0
      ccsize = 1
   100 lbegin = lvlend + 1
      lvlend = ccsize
      do 400 i = lbegin, lvlend
         node = ls(i)
         jstrt = -xadj(node)
         jstop = iabs(xadj(node + 1)) - 1
         ideg = 0
         if (jstop.lt.jstrt) go to 300
         do 200 j = jstrt, jstop
            nbr = adjncy(j)
            if (mask(nbr).eq.0) go to 200
            ideg = ideg + 1
            if (xadj(nbr).lt.0) go to 200
            xadj(nbr) = -xadj(nbr)
            ccsize = ccsize + 1
            ls(ccsize) = nbr
   200   continue
   300   deg(node) = ideg
   400 continue
      lvsize = ccsize - lvlend
      if (lvsize.gt.0) go to 100
      do 500 i = 1, ccsize
         node = ls(i)
         xadj(node) = -xadj(node)
   500 continue
   
   
   end

   
!c###############################################!c

!      subroutine type_arete(i, i1, i2, ne_d1, nu_d1, typ_el_d1)

!     use numbatmod
!     integer i, i1, i2, ne_d1
!     integer nu_d1(3,ne_d1), typ_el_d1(ne_d1)
!     integer j, k, k1, k2

!     i = 0
!     do j=1,ne_d1
!       k = typ_el_d1(j)
!       k1 = nu_d1(1,j)
!       k2 = nu_d1(2,j)
!       if(k1 .eq. i1 .and. k2 .eq. i2) then
!         i = k
!         return
!       elseif(k1 .eq. i2 .and. k2 .eq. i1) then
!         i = k
!         return
!       endif
!     enddo
!     print*, '?? type_arete: Resultat negatif'
!     print*, 'i1, i2 = ', i1, i2
!     stop

!     end

!    ------------------------------------------------------------------
!        subroutine mailp2(tcp2,maxvis,ne,n_pts,ns,liste,nb,numero)

! !    ------------------------------------------------------------------

!       integer flag, i, iplusj_mod3(3), is1, is2, j, jel, jj, ne, n_pts, ns, temp
!       integer liste(maxvis,ns), maxvis, nb(ns), numero(maxvis,ns)
!       integer tcp2(6,ne)
! c
! !    ------------------------------------------------------------------
! c
! !     print*, 'MAILP2: 0: n_pts = ', n_pts
!       iplusj_mod3(1) = 2
!       iplusj_mod3(2) = 3
!       iplusj_mod3(3) = 1
! c
!       do i = 1, ns
!          nb(i) = 0
!       enddo
! c
!       n_pts = ns
!       do jel = 1, ne
!          do i = 1, 3
! c
!             is1 = tcp2(i,jel)
!             is2 = tcp2(iplusj_mod3(i),jel)
!             if (is1.gt.is2) then
!                temp = is1
!                is1  = is2
!                is2  = temp
!             endif
! c
!             flag = 0
! c
!             if (nb(is1).eq.0) go to 1
! c
!             do j = 1, nb(is1)
!                if (liste(j,is1).eq.is2) then
!                   flag = 1
!                   jj   = j
!                   goto 1
!                endif
!             enddo
! c
! 1           continue
!             if (flag.eq.0) then
! c
! !             l'arete n'est pas dans la liste
! !             -------------------------------
! c
!                n_pts = n_pts+1
!                tcp2(i+3,jel) = n_pts
!                nb(is1) = nb(is1)+1
!                liste(nb(is1),is1) = is2
!                numero(nb(is1),is1) = n_pts
!             else
! c
! !             l'arete est deja dans la liste
! !             ------------------------------
! c
!                tcp2(i+3,jel) = numero(jj,is1)
!             endif
! c
!          enddo
! c
!       enddo
! c
!       end





!    ------------------------------------------------------------------!!     
! subroutine prepmailp2(tcp1,maxvis,ne,ns,visited)!c!    ------------------------------------------------------------------!
!     integer i, ii, jel, maxvis, ne, ns, tcp1(6,ne), visited(ns)!c!    ------------------------------------------------------------------!c
!     do i = 1, ns
!        visited(i) = 0
!     enddo!!     do jel = 1, ne
!        do i = 1, 3
!           ii = tcp1(i,jel)
!           visited(ii) = visited(ii)+1
!        enddo
!     enddo!!     maxvis = 0
!     do i = 1, ns
!        maxvis = max(maxvis,visited(i))
!     enddo!c
!     end!cc
!    ------------------------------------------------------------------



!*********************************************************************


!     subroutine test_border(n_pts, type_nod0, x, y, i_err,
!    *                       file1_mesh, ui)

!     use numbatmod
!     integer n_pts, i_err
!     integer type_nod0(n_pts)
!     double precision x(n_pts), y(n_pts)

!     integer type_nod, n_typ2
!     double precision d_period, tmp, tmp2

!     integer  n_border1, n_border2, n_typ1, ui
!     integer i, j, i1, i2, j1, j2, max_period
!     parameter(max_period=2000)
!     double precision x2(2), x1(2), y_min, y_max
!     double precision delta_x, tol
!     integer period1(max_period), period2(max_period)
!     integer period3(2,max_period)
!     character*(*) file1_mesh
!     character file_ui*100
!     common/imp_file/file_ui

!     character fichier_mail*20

!    debut de la programmation de la sous-routine maillage
!    -----------------------------------------------------

!     i_err = 0
!     tol = 1.0d-6
!     y_min = y(1)
!     y_max = y(1)
!     do i=1,n_pts
!       if(y(i) .lt. y_min) y_min = y(i)
!       if(y(i) .gt. y_max) y_max = y(i)
!     enddo
!     i1 = 0
!     i2 = 0
!     do i=1,n_pts
!       if(abs(y(i)-y_max) .lt. tol) then
!         if(i1 .eq. 0) then
!           i1 = 1
!           x1(1) = x(i)
!           x1(2) = x(i)
!         endif
!         if(x(i) .lt. x1(1)) x1(1) = x(i)
!         if(x(i) .gt. x1(2)) x1(2) = x(i)
!       endif
!       if(abs(y(i)-y_min) .lt. tol) then
!         if(i2 .eq. 0) then
!           i2 = 1
!           x2(1) = x(i)
!           x2(2) = x(i)
!         endif
!         if(x(i) .lt. x2(1)) x2(1) = x(i)
!         if(x(i) .gt. x2(2)) x2(2) = x(i)
!       endif
!     enddo
!     delta_x = x1(1) - x2(1)
!     d_period = x1(2)-x1(1)

!     n_border1 = 0
!     n_border2 = 0
!     n_typ1 = 0
!     n_typ2 = 0

!     do i=1,n_pts
!       type_nod = type_nod0(i)
!       if(type_nod .eq. 1) then
!         n_typ1 = n_typ1 + 1
!       elseif(type_nod .eq. 2) then
!         n_typ2 = n_typ2 + 1
!       elseif(type_nod .eq. 3) then
!         n_border1 = n_border1 + 1
!        period1(n_border1) = i
!       elseif(type_nod .eq. 4) then
!         n_border2 = n_border2 + 1
!        period2(n_border2) = i
!       endif
!      if(n_border2 .ge. max_period) then
!       open (unit=ui,file=file_ui)
!          write(*,*) 'get_border (in conv_gmsh.f): attention'
!          write(*,*) 'n_border2 >= max_period = ',
!    *        n_border2, max_period
!       close(ui)
!        i_err = 1
!        stop
!      endif
!     enddo

!     if(n_border1 .ne. n_border2) then
!       open (unit=ui,file=file_ui)
!         write(*,*) 'get_border (in conv_gmsh.f): ',
!    *    ' n_border1 .ne. n_border2'
!       close(ui)
!         i_err = 2
!       stop
!     endif
!     do i=1,n_border1
!       i1 = period1(i)
!       period3(1,i) = i1
!       period3(2,i) = 0
!       do j=1,n_border1
!         j1=period2(j)
!         tmp=dsqrt((x(j1)-x(i1))**2+(y(j1)-y(i1))**2)
!         tmp2 = abs(y(j1)-y(i1))
!       if(abs(tmp-d_period) .lt. 5.0d-6 .and. tmp2 .lt. 5.0d-6) then
!           if(period3(2,i) .ne. 0) then
!          open (unit=ui,file=file_ui)
!             write(*,*) 'get_border  (in conv_gmsh.f): ',
!    *        'i,j = ', i, j, i1, j1
!             write(*,*) 'probleme ave!period3(2,i) = ', period3(2,i)
!             write(*,*) x(i1), y(i1), d_period
!             write(*,*) x(j1), y(j1), tmp, (tmp-d_period)
!             j2 = period3(2,i)
!             tmp=dsqrt((x(j2)-x(i1))**2+(y(j2)-y(i1))**2)
!             write(*,*) x(j2), y(j2), tmp, (tmp-d_period)
!             i_err = 3
!             close(ui)
!             stop
!           endif
!           period3(2,i) = j1
!         endif
!       enddo
!     enddo
!     do i=1,n_border1
!       if(period3(2,i) .eq. 0) then
!          open (unit=ui,file=file_ui)
!             write(*,*) 'get_border (in conv_gmsh.f): ',
!    *        ' period3(2,i) = 0', i
!          close(ui)
!         i_err = 4
!         stop
!       endif
!     enddo

!     
!     end

!    ------------------------------------------------------------------
