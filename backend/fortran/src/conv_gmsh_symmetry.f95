!c###############################################!!c###############################################
subroutine symmetry(n_pts, ne, max_n_gelts_triangs, idfn, nu, typ_el, x, y, i_sym)

   !!    symmetry: symmetrize an FEM mesh!c*******************************************************

   use numbatmod
   integer i_sym, max_n_gelts_triangs

   integer ne, n_pts
   integer nu(6,max_n_gelts_triangs), typ_el(max_n_gelts_triangs)
   integer idfn(MAX_N_PTS)
   double precision x(MAX_N_PTS), y(MAX_N_PTS)

   integer max_n_gelts_triangs_0, MAX_N_PTS_0
   parameter(MAX_N_PTS_0=250000)
   parameter (max_n_gelts_triangs_0=120000)
   integer ne_0, n_pts_0, idfn_0(MAX_N_PTS_0)
   integer nu_0(6,max_n_gelts_triangs_0), typ_el_0(max_n_gelts_triangs_0)
   double precision x_0(MAX_N_PTS_0),  y_0(MAX_N_PTS_0)
   integer tab_ne(max_n_gelts_triangs_0), tab_n_pts(MAX_N_PTS_0,3)!c
   integer i, j

   n_pts_0 = n_pts
   ne_0 = ne
   do i=1,n_pts_0
      x_0(i) = x(i)
      y_0(i) = y(i)
      idfn_0(i) = idfn(i)
   enddo
   do i=1,ne_0
      do j=1,6
         nu_0(j,i) = nu(j,i)
      enddo
      typ_el_0(i) = typ_el(i)
   enddo!ccccccccccc
   if(i_sym .eq. 1) then
      call y_symmetry(n_pts, ne, ne_0, n_pts_0, max_n_gelts_triangs, idfn, nu, typ_el, &
         idfn_0, nu_0, typ_el_0, x, y, x_0, y_0, tab_ne, tab_n_pts)

   elseif(i_sym .eq. 2) then
      call x_symmetry(n_pts, ne, ne_0, n_pts_0, max_n_gelts_triangs, idfn, nu, typ_el, &
         idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,    tab_ne, tab_n_pts)
   elseif(i_sym .eq. 3) then
      call y_symmetry(n_pts, ne, ne_0, n_pts_0, max_n_gelts_triangs, idfn, nu, typ_el, &
         idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,  tab_ne, tab_n_pts)

      n_pts_0 = n_pts
      ne_0 = ne
      do i=1,n_pts_0
         x_0(i) = x(i)
         y_0(i) = y(i)
         idfn_0(i) = idfn(i)
      enddo
      do i=1,ne_0
         do j=1,6
            nu_0(j,i) = nu(j,i)
         enddo
         typ_el_0(i) = typ_el(i)
      enddo

      call x_symmetry(n_pts, ne, ne_0, n_pts_0, max_n_gelts_triangs, idfn, nu, typ_el, idfn_0, &
         nu_0, typ_el_0, x, y, x_0, y_0, tab_ne, tab_n_pts)

   else
      return
   endif



end

 !c###############################################

subroutine y_symmetry(n_pts, ne, ne_0, n_pts_0,  max_n_gelts_triangs, idfn, nu, typ_el, &
   idfn_0, nu_0, typ_el_0, x, y, x_0, y_0, tab_ne, tab_n_pts)
   use numbatmod
   integer max_n_gelts_triangs

   integer ne_0, n_pts_0, idfn_0(MAX_N_PTS)
   integer nu_0(6,max_n_gelts_triangs), typ_el_0(max_n_gelts_triangs)

   integer ne, n_pts, idfn(MAX_N_PTS)
   integer nu(6,max_n_gelts_triangs), typ_el(max_n_gelts_triangs)

   integer tab_ne(max_n_gelts_triangs), tab_n_pts(MAX_N_PTS,3)
   double precision x(MAX_N_PTS), y(MAX_N_PTS)
   double precision x_0(MAX_N_PTS),  y_0(MAX_N_PTS)!!    Local variables
   integer i, i1, i2, i_a, i_b, j, j1, j2
   integer ne_1, n_pts_1, n_pts_2
   integer t_nodes_0(6), t_nodes_a(6), t_nodes_b(6)
   double precision tol, y_min, y_max, y_mid
   double precision x_a, y_a, x_b, y_b

   integer ui, debug
   character file_ui*100
   common/imp/ui, debug
   common/imp_file/file_ui
   y_min = y_0(1)
   y_max = y_0(1)
   do i=1,n_pts_0
      if(y_0(i) .lt. y_min) y_min = y_0(i)
      if(y_0(i) .gt. y_max) y_max = y_0(i)
   enddo

   !    Selecting points in upper half of the domain
   tol = 1.0d-7
   y_mid = (y_min+y_max)/2.0d0
   i1 = 0
   i2 = 0
   do i=1,n_pts_0
      if(y_0(i) .ge. y_mid-tol) then
         i1 = i1 + 1
         tab_n_pts(i,1) = i1
         if(Abs(y_0(i)-y_mid) .le. tol) then
            !          No duplication for a point on the symmetry line
            i2 = i2 + 1
            tab_n_pts(i,2) = 1
         else
            tab_n_pts(i,2) = 2
         endif
      else
         tab_n_pts(i,1) = 0
         tab_n_pts(i,2) = 0
      endif
   enddo
   n_pts_1 = i1
   n_pts_2 = i2
   n_pts = 2*i1 - i2!!    Selecting triangles in upper half of the domain
   i1 = 0
   do i=1,ne_0
      i2 = 0
      do j=1,3
         j1 = nu_0(j,i)
         j2 = tab_n_pts(j1,1)
         if(j2 .gt. 0) i2 = i2 + 1
      enddo
      if(i2 .eq. 3) then
         i1 = i1 + 1
         tab_ne(i) = i1
      else
         tab_ne(i) = 0
      endif
   enddo
   ne_1 = i1
   ne = 2*i1
   !!    Generating the symmetrized FEM mesh
   i_b = n_pts_1
   do i=n_pts_0,1,-1
      tab_n_pts(i,3) = 0
      x_a = x_0(i)
      y_a = y_0(i)

      x_b = x_a
      y_b = 2.0d0*y_mid - y_a

      i_a = tab_n_pts(i,1)
      !       i_b = n_pts - i_a + 1
      i1 = tab_n_pts(i,2)

      if(i_a .gt. 0) then
         x(i_a) = x_a
         y(i_a) = y_a
         idfn(i_a) = idfn_0(i)
         if(i1 .eq. 2) then
            i_b = i_b + 1
            tab_n_pts(i,3) = i_b
            x(i_b) = x_b
            y(i_b) = y_b
            idfn(i_b) = idfn_0(i)
            if(idfn(i_b) .eq. 1) idfn(i_b) = 2
         endif
      endif
   enddo!ccccccccccccccccccccccc
   i_b = ne_1
   do i=ne_0,1,-1
      i_a = tab_ne(i)
      if(i_a .gt. 0) then
         i_b = i_b + 1
         do j=1,6
            t_nodes_0(j) = nu_0(j,i)
         enddo

         do j=1,6
            j1 = tab_n_pts(t_nodes_0(j),1)
            t_nodes_a(j) = j1
            j2 = tab_n_pts(t_nodes_0(j),2)

            if(j2 .eq. 1) then
               t_nodes_b(j) = j1
            elseif(j2 .eq. 2) then
               t_nodes_b(j) = tab_n_pts(t_nodes_0(j),3)
               !             t_nodes_b(j) = n_pts - j1 + 1
            else
               open (unit=ui,file=file_ui)
               write(*,*) 'SYMMETRY: tab_n_pts(i,2) = ', j2
               write(*,*) 'i, tab_n_pts(i,1) = ', t_nodes_0(j), j1
               stop
            endif
         enddo
         do j=1,6
            nu(j,i_a) = t_nodes_a(j)
         enddo
         typ_el(i_a) = typ_el_0(i)

         !         i_b = ne - i_a + 1
         if(i_a .gt. ne/2) stop
         do j=1,3
            nu(j,i_b) = t_nodes_b(3-j+1)
            nu(j+3,i_b) = t_nodes_b(6-j+1)
         enddo!!      Symmetry reverses the orientation
         !      so we must reverse the numbering to get the positive orientation

         nu(1,i_b) = t_nodes_b(1)
         nu(2,i_b) = t_nodes_b(3)
         nu(3,i_b) = t_nodes_b(2)
         nu(4,i_b) = t_nodes_b(6)
         nu(5,i_b) = t_nodes_b(5)
         nu(6,i_b) = t_nodes_b(4)

         typ_el(i_b) = typ_el_0(i)
      endif
   enddo
end

 !c###############################################!!c###############################################

subroutine x_symmetry(n_pts, ne, ne_0, n_pts_0, max_n_gelts_triangs, idfn, nu, typ_el, &
   idfn_0, nu_0, typ_el_0, x, y, x_0, y_0, tab_ne, tab_n_pts)
   use numbatmod
   integer max_n_gelts_triangs

   integer ne_0, n_pts_0, idfn_0(MAX_N_PTS)
   integer nu_0(6,max_n_gelts_triangs), typ_el_0(max_n_gelts_triangs)

   integer ne, n_pts, idfn(MAX_N_PTS)
   integer nu(6,max_n_gelts_triangs), typ_el(max_n_gelts_triangs)

   integer tab_ne(max_n_gelts_triangs), tab_n_pts(MAX_N_PTS,3)
   double precision x(MAX_N_PTS), y(MAX_N_PTS)
   double precision x_0(MAX_N_PTS),  y_0(MAX_N_PTS)!!    Local variables
   integer i, i1, i2, i_a, i_b, j, j1, j2
   integer ne_1, n_pts_1, n_pts_2
   integer t_nodes_0(6), t_nodes_a(6), t_nodes_b(6)
   double precision tol, x_min, x_max, x_mid
   double precision x_a, y_a, x_b, y_b

   integer ui, debug
   character file_ui*100
   common/imp/ui, debug
   common/imp_file/file_ui
   x_min = x_0(1)
   x_max = x_0(1)
   do i=1,n_pts_0
      if(x_0(i) .lt. x_min) x_min = x_0(i)
      if(x_0(i) .gt. x_max) x_max = x_0(i)
   enddo

   !    Selecting points in upper half of the domain
   tol = 1.0d-7
   x_mid = (x_min+x_max)/2.0d0
   i1 = 0
   i2 = 0
   do i=1,n_pts_0
      if(x_0(i) .le. x_mid+tol) then
         i1 = i1 + 1
         tab_n_pts(i,1) = i1
         if(Abs(x_0(i)-x_mid) .le. tol) then
            !          No duplication for a point on the symmetry line
            i2 = i2 + 1
            tab_n_pts(i,2) = 1
         else
            tab_n_pts(i,2) = 2
         endif
      else
         tab_n_pts(i,1) = 0
         tab_n_pts(i,2) = 0
      endif
   enddo
   n_pts_1 = i1
   n_pts_2 = i2
   n_pts = 2*i1 - i2!!    Selecting triangles in left half of the domain
   i1 = 0
   do i=1,ne_0
      i2 = 0
      do j=1,3
         j1 = nu_0(j,i)
         j2 = tab_n_pts(j1,1)
         if(j2 .gt. 0) i2 = i2 + 1
      enddo
      if(i2 .eq. 3) then
         i1 = i1 + 1
         tab_ne(i) = i1
      else
         tab_ne(i) = 0
      endif
   enddo
   ne_1 = i1
   ne = 2*i1
   !!    Generating the symmetrized FEM mesh
   i_b = n_pts_1
   do i=n_pts_0,1,-1
      tab_n_pts(i,3) = 0
      x_a = x_0(i)
      y_a = y_0(i)

      x_b = 2.0d0*x_mid - x_a
      y_b = y_a

      i_a = tab_n_pts(i,1)
      i1 = tab_n_pts(i,2)

      if(i_a .gt. 0) then
         x(i_a) = x_a
         y(i_a) = y_a
         idfn(i_a) = idfn_0(i)
         if(i1 .eq. 2) then
            i_b = i_b + 1
            tab_n_pts(i,3) = i_b
            x(i_b) = x_b
            y(i_b) = y_b
            idfn(i_b) = idfn_0(i)
            if(idfn(i_b) .eq. 3) idfn(i_b) = 4
         endif
      endif
   enddo!ccccccccccccccccccccccc!         i_b = ne_1
   do i=ne_0,1,-1
      i_a = tab_ne(i)
      if(i_a .gt. 0) then
         i_b = i_b + 1
         do j=1,6
            t_nodes_0(j) = nu_0(j,i)
         enddo

         do j=1,6
            j1 = tab_n_pts(t_nodes_0(j),1)
            t_nodes_a(j) = j1
            j2 = tab_n_pts(t_nodes_0(j),2)

            if(j2 .eq. 1) then
               t_nodes_b(j) = j1
            elseif(j2 .eq. 2) then
               t_nodes_b(j) = tab_n_pts(t_nodes_0(j),3)
               !             t_nodes_b(j) = n_pts - j1 + 1
            else
               open (unit=ui,file=file_ui)
               write(*,*) 'SYMMETRY_X: tab_n_pts(i,2) = ', j2
               write(*,*) 'i, tab_n_pts(i,1) = ', t_nodes_0(j), j1
               stop
            endif
         enddo
         do j=1,6
            nu(j,i_a) = t_nodes_a(j)
         enddo
         typ_el(i_a) = typ_el_0(i)

         !         i_b = ne - i_a + 1
         if(i_a .gt. ne/2) stop
         do j=1,3
            nu(j,i_b) = t_nodes_b(3-j+1)
            nu(j+3,i_b) = t_nodes_b(6-j+1)
         enddo!!      Symmetry reverses the orientation
         !      so we must reverse the numbering to get the positive orientation

         nu(1,i_b) = t_nodes_b(1)
         nu(2,i_b) = t_nodes_b(3)
         nu(3,i_b) = t_nodes_b(2)
         nu(4,i_b) = t_nodes_b(6)
         nu(5,i_b) = t_nodes_b(5)
         nu(6,i_b) = t_nodes_b(4)

         typ_el(i_b) = typ_el_0(i)
      endif
   enddo
end

