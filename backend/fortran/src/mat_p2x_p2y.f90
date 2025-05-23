!
!!!!!!!!!!!!!!!!
!
!      Compute the integral: m_ij = Integrate[D[P(i),x] * D[P(j),y], over a triangle]
!
!!!!!!!!!!!!!!!!
!
subroutine find_overlaps_p2x_p2y (mat, b_inv_t, det_b)
!
!!!!!!!!!!!!!!!!
!
   implicit none
   double precision mat(6,6), b_inv_t(2,2), det_b

   double precision dt, factor, factor2
   double precision c11, c12, c21, c22
   integer(8) i, j

   dt=dabs(det_b)
   c11 = b_inv_t(1, 1)
   c12 = b_inv_t(1, 2)
   c21 = b_inv_t(2, 1)
   c22 = b_inv_t(2, 2)

   mat(1,1) = 3*(c11 + c12)*(c21 + c22)
   mat(1,2) = (c11 + c12)*c21
   mat(1,3) = (c11 + c12)*c22
   mat(1,4) = -4*(c11 + c12)*c21
   mat(1,5) = 0
   mat(1,6) = -4*(c11 + c12)*c22
   mat(2,1) = c11*(c21 + c22)
   mat(2,2) = 3*c11*c21
   mat(2,3) = -(c11*c22)
   mat(2,4) = -4*c11*(c21 + c22)
   mat(2,5) = 4*c11*c22
   mat(2,6) = 0
   mat(3,1) = c12*(c21 + c22)
   mat(3,2) = -(c12*c21)
   mat(3,3) = 3*c12*c22
   mat(3,4) = 0
   mat(3,5) = 4*c12*c21
   mat(3,6) = -4*c12*(c21 + c22)
   mat(4,1) = -4*c11*(c21 + c22)
   mat(4,2) = -4*(c11 + c12)*c21
   mat(4,3) = 0
   mat(4,4) = 4*(c11*(2*c21 + c22) + c12*(c21 + 2*c22))
   mat(4,5) = -4*(c11*c22 + c12*(c21 + 2*c22))
   mat(4,6) = 4*(c12*c21 + c11*c22)
   mat(5,1) = 0
   mat(5,2) = 4*c12*c21
   mat(5,3) = 4*c11*c22
   mat(5,4) = -4*(c11*c22 + c12*(c21 + 2*c22))
   mat(5,5) = 4*(c11*(2*c21 + c22) + c12*(c21 + 2*c22))
   mat(5,6) = -4*(c12*c21 + c11*(2*c21 + c22))
   mat(6,1) = -4*c12*(c21 + c22)
   mat(6,2) = 0
   mat(6,3) = -4*(c11 + c12)*c22
   mat(6,4) = 4*(c12*c21 + c11*c22)
   mat(6,5) = -4*(c12*c21 + c11*(2*c21 + c22))
   mat(6,6) = 4*(c11*(2*c21 + c22) + c12*(c21 + 2*c22))


   factor = 6
   factor2 = dt / factor
   do j=1,6
      do i=1,6
         mat(i,j) = mat(i,j) * factor2
      enddo
   enddo
!
   return
end

