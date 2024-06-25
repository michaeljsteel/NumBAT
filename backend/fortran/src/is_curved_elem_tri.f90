
!
!     Check if an element (triangle) has a curved face
!
!  Currently always says no!




subroutine is_curved_elem_tri_impl (nnodes, xel, info_curved, tmp)
!
   implicit none
   integer(8) nnodes, info_curved
   double precision xel(2,nnodes)
   double precision tmp


   integer(8) nnd_triangle
   parameter (nnd_triangle = 6)
   double precision xel_triangle(2,nnd_triangle)

   integer(8) i, j, i2
   double precision tol

   info_curved = 0

   if (nnodes .ne. nnd_triangle) then
      write(*,*)
      write(*,*) "   ???"
      write(*,*) "is_curved_elem_tri: nnodes != nnd_triangle : ",&
      &nnodes, nnd_triangle
      write(*,*) "is_curved_elem_tri: Aborting..."
      stop
   endif

!     Vertices
   do i=1,3
      do j=1,2
         xel_triangle(j,i) = xel(j,i)
      enddo
   enddo

!     Mid-points
   do i=1,3
      i2 = modulo(i+1,3_8)
      if(i2 .eq. 0) i2 = 3
      do j=1,2
         xel_triangle(j,i+3) = (xel(j,i)+xel(j,i2))/2.0d0
      enddo
   enddo

   tmp = 0.0d0
   do i=1,nnodes
      do j=1,2
         tmp = tmp + (xel_triangle(j,i) - xel(j,i))**2
      enddo
   enddo

   tol = 1.0d-14
   if(abs(tmp) .lt. tol) then
      info_curved = 0
   else
      info_curved = 1
   endif


   return
end
