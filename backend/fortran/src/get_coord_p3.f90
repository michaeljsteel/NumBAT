!!!!!!!!!!!!!!!!

!  Input : initial P2 FEM mesh
!  Output: Set the coordinates and node type of the P3 node

!  Fills x_N_E_F and type_N_E_F arrays
!!!!!!!!!!!!!!!!

subroutine get_coord_p3(n_msh_el, n_msh_pts, nodes_per_el, n_ddl, &
   elnd_to_mesh, type_nod, table_N_E_F, &
   type_N_E_F, v_nd_xy, x_N_E_F, visited)


   implicit none
   integer(8) n_msh_el, n_msh_pts, nodes_per_el, n_ddl
   integer(8) elnd_to_mesh(nodes_per_el,n_msh_el), table_N_E_F(14,n_msh_el)
   integer(8) type_nod(n_msh_pts), type_N_E_F(2,n_ddl)
   integer(8) visited(n_ddl)

   double precision v_nd_xy(2,n_msh_pts), x_N_E_F(2,n_ddl)

   integer(8) nddl_0
   parameter (nddl_0 = 14)
   integer(8) j, k, k1, n, ind, ip(2,3)
   integer(8) iel, inod, inod1, inod2, mm
   integer(8) nut0(6), nut_N_E_F(nddl_0)

   double precision xx1, xx2, xx3, yy1, yy2, yy3
   double precision dx1, dy1
   double precision tmp1, tmp2, tmp3, tmp_x, tmp_y


   if ( nodes_per_el .ne. 6 ) then
      write(*,*) "get_coord_p3: problem nodes_per_el = ", nodes_per_el
      write(*,*) "get_coord_p3: nodes_per_el should be equal to 6 !"
      write(*,*) "get_coord_p3: Aborting..."
      stop
   endif

!  ip(1,i) = i+1 MOD 3 ; Number of the next vertex to vertex i
!  ip(2,i) = i+2 MOD 3 ; Number of the second next vertex to vertex i
   ip(1,1) = 2
   ip(1,2) = 3
   ip(1,3) = 1

   ip(2,1) = 3
   ip(2,2) = 1
   ip(2,3) = 2

   do j=1,n_ddl
      visited(j) = 0
   enddo


   !  The first 4 entries of table_N_E_F(*,i) correspond to face and edges and have been done
   mm = 4
   do iel=1,n_msh_el

      do inod=1,nodes_per_el
         nut0(inod) = elnd_to_mesh(inod,iel)
      enddo

      !  the 10 node of a P3 element
      do inod=1,10
         nut_N_E_F(inod) = table_N_E_F(inod+mm,iel)
      enddo

      !  scan the vertices ############
      do inod=1,3
         k = nut0(inod)
         if(visited(k) .eq. 0) then
            visited(k) = iel
            inod1 = nut0(inod)
            inod2 = nut_N_E_F(inod)
            x_N_E_F(1,inod2) = v_nd_xy(1,inod1)
            x_N_E_F(2,inod2) = v_nd_xy(2,inod1)
            type_N_E_F(1,inod2) = type_nod(inod1)

            !  Vertex => dimension zero
            type_N_E_F(2,inod2) = 0
         endif
      enddo

      !  scan the nodes located on the edges ############
      do inod=4,nodes_per_el

         k=nut0(inod)
         if(k .lt. 1) then
            print*, 'k = ', k
            print*, 'iel, inod = ', iel, inod
            print*, 'nut0 = ', (nut0(inod2), inod2=1,nodes_per_el)
            stop
         endif

         if(visited(k) .eq. 0) then
            visited(k) = iel
!  Endpoints of the edge
            k1 = nut0(inod-3)
            xx1 = v_nd_xy(1,k1)
            yy1 = v_nd_xy(2,k1)
            k1 = nut0(ip(1,inod-3))
            xx2 = v_nd_xy(1,k1)
            yy2 = v_nd_xy(2,k1)
            dx1 = (xx2-xx1)/3.0d0
            dy1 = (yy2-yy1)/3.0d0

            !  type of the mid-edge node of the initial P2 mesh
            ind = type_nod(nut0(inod))

            !  2 nodes per edge (for P3 element)
            do inod2=1,2
               k1 = nut_N_E_F(inod2+2*(inod-4)+3)
               x_N_E_F(1,k1) = xx1 + inod2*dx1
               x_N_E_F(2,k1) = yy1 + inod2*dy1
               type_N_E_F(1,k1) = ind

               !  Node => dimension zero
               type_N_E_F(2,k1) = 0
            enddo
         endif
      enddo

      !  Coordinate of the vertices
      k1 = nut0(1)
      xx1 = v_nd_xy(1,k1)
      yy1 = v_nd_xy(2,k1)
      k1 = nut0(2)
      xx2 = v_nd_xy(1,k1)
      yy2 = v_nd_xy(2,k1)
      k1 = nut0(3)
      xx3 = v_nd_xy(1,k1)
      yy3 = v_nd_xy(2,k1)

      !  The tenth node is at the center of the triangle
      !  dimension(P3) = 10
      n = 10

      !  this node is an interior node of the triangle ############
      k1 = nut_N_E_F(n)
      tmp1 = 1.0d0/3.0d0
      tmp2 = 1.0d0/3.0d0
      tmp3 = 1.0d0/3.0d0
      tmp_x = xx1*tmp1+xx2*tmp2+xx3*tmp3  !  TODO: clean me
      tmp_y = yy1*tmp1+yy2*tmp2+yy3*tmp3
      x_N_E_F(1,k1) = tmp_x
      x_N_E_F(2,k1) = tmp_y

      !  interior node
      type_N_E_F(1,k1) = 0

      !  Node => dimension zero
      type_N_E_F(2,k1) = 0

   enddo

   return
end
