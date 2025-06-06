 !
 !***********************************************************************
 !
subroutine periodic_node (nel, npt, nnodes, type_nod, x,&
&ip_period, n_period, m_elnd_to_mshpt, lat_vecs)
   !
   !***********************************************************************
   !
   !   type_nod = 0  => interiour point
   !   type_nod != 0 => boundary point
   !
   !***********************************************************************
   !
   implicit none
   integer(8) nel, npt, nnodes
   integer(8) type_nod(npt), ip_period(npt)
   integer(8) n_period(npt)
   integer(8) m_elnd_to_mshpt(nnodes,nel)
   double precision lat_vecs(2,2)
   double precision x(2,npt)
   integer(8) i, j, i1, j1, i_not_periodic
   double precision period_x, period_y
   double precision x_min, y_min
   double precision x_max, y_max
   double precision tmp1, tmp2, tol
   double precision x_r, y_r
   double precision delta_v(2),  vec(2)
   integer(8) ix, iy, test_lattice
   integer(8) list_end(2,3), j2, k, debug
   !
   debug = 0
   x_min = x(1,1)
   x_max = x(1,1)
   do i=1,npt
      x_r = x(1,i)
      if(x_r .lt. x_min) x_min = x_r
      if(x_r .gt. x_max) x_max = x_r
   enddo
   y_min = x(2,1)
   y_max = x(2,1)
   do i=1,npt
      y_r = x(2,i)
      if(y_r .lt. y_min) y_min = y_r
      if(y_r .gt. y_max) y_max = y_r
   enddo
   if (debug .eq. 1) then
      write(*,*) "periodic_node: x_min, x_max = ", x_min, x_max
      write(*,*) "periodic_node: y_min, y_max = ", y_min, y_max
   endif
   !
   tol = 1.0d-6
   period_x = x_max - x_min
   period_y = y_max - y_min
   do i=1,npt
      ip_period(i) = 0
      n_period(i) = 0
   enddo
   do i=1,npt-1
      if(type_nod(i) .ne. 0) then
         do j=i+1,npt
            if(type_nod(j) .ne. 0  .and. ip_period(j) .eq. 0) then
               delta_v(1) = x(1,j) - x(1,i)
               delta_v(2) = x(2,j) - x(2,i)
               test_lattice = 0
               do ix=-1,1
                  do iy=-1,1
                     do k=1,2
                        vec(k) = ix*lat_vecs(k,1)&
                        &+ iy*lat_vecs(k,2)
                     enddo
                     tmp1 = 0.0d0
                     do k=1,2
                        tmp2 = delta_v(k) - vec(k)
                        tmp1 = tmp1 + abs(tmp2)
                     enddo
                     if (tmp1 .lt. tol) then
                        test_lattice = 1
                        goto 10
                     endif
                     !                  enddo
                  enddo
               enddo
10             continue
               !
               if(test_lattice .eq. 1) then
                  n_period(i) = n_period(i) + 1
                  i1 = ip_period(i)
                  j1 = ip_period(j)
                  if(i1 .eq. 0) then
                     ip_period(i) = i
                     ip_period(j) = i
                  elseif(i1 .ne. 0 .and. j1 .eq. 0) then
                     ip_period(j) = i1
                  else
                     write(*,*)
                     write(*,*) "  ???"
                     write(*,*) "periodic_node: for i, j = ", i, j
                     write(*,*) "periodic_node: ip_period : ", i1, j1
                     write(*,*) "periodic_node: Aborting..."
                     stop
                  endif
               endif
            endif
         enddo
      endif
   enddo
   !
   i_not_periodic = 0
   do i=1,npt
      if(type_nod(i) .ne. 0 .and. ip_period(i) .eq. 0) then
         if(i_not_periodic .eq. 0) then
            open (unit=11,file="not_period.txt",status='unknown')
         endif
         i_not_periodic = i_not_periodic + 1
         write(11,*) i, type_nod(i), ip_period(i), n_period(i)
      endif
   enddo

   if(i_not_periodic .gt. 0) then
      close ( unit = 11)
      write(*,*)
      write(*,*) "  ???"
      write(*,*) "periodic_node: the FEM mesh is not periodic"
      write(*,*) "periodic_node: see the file not_period.txt"
      write(*,*) "periodic_node: Aborting..."
      write(*,*)
      stop
   endif
   !
   if (debug .eq. 1) then
      open (unit = 10, file="ip_period.txt",status='unknown')
      do i=1,npt
         write(10,*) i, type_nod(i), ip_period(i), n_period(i),&
         &"     ", (x(j,i),j=1,2)
      enddo
      close ( unit = 10)
   endif
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !     some check ...
   !
   !     Endpoints of the 6 edges (mid-point) of the reference triangle
   !
   i = 1
   list_end(1,i) = 1
   list_end(2,i) = 2
   i = 2
   list_end(1,i) = 2
   list_end(2,i) = 3
   i = 3
   list_end(1,i) = 3
   list_end(2,i) = 1


   do i=1,nel
      do j=4,nnodes
         j1 = list_end(1,j-3)
         j2 = list_end(2,j-3)
         if(ip_period(m_elnd_to_mshpt(j,i)) .ne. 0) then
            if (ip_period(m_elnd_to_mshpt(j1,i)) .eq. 0 .or.&
            &ip_period(m_elnd_to_mshpt(j2,i)) .eq. 0) then
               write(*,*) "periodic_node: m_elnd_to_mshpt = ",&
               &type_nod(m_elnd_to_mshpt(j1,i)), type_nod(m_elnd_to_mshpt(j2,i)),&
               &type_nod(m_elnd_to_mshpt(j,i))
               write(*,*) "periodic_node: m_elnd_to_mshpt(k,i): ",&
               &(m_elnd_to_mshpt(k,i),k=1,6)
               write(*,*) "periodic_node: ip_period(m_elnd_to_mshpt(k,i)): ",&
               &(ip_period(m_elnd_to_mshpt(k,i)),k=1,6)
               write(*,*) "periodic_node: type_nod(m_elnd_to_mshpt(k,i)): ",&
               &(type_nod(m_elnd_to_mshpt(k,i)),k=1,6)
               write(*,*) "periodic_node: Aborting..."
               stop
            endif
         endif
      enddo
   enddo

   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   return
end
 !
