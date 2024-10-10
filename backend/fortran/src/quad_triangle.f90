
subroutine quad_triangle (nquad, nquad_max, wq, xq, yq)

   !***********************************************************************
   !
   ! Evaluates the elementary integrals of the convective components on each triangle.
   ! Here, the 16-point Gauss-Hammer method is used,which exactly integrates polynomials of the eighth degree.
   !
   !       Reference
   !       J. N. Lyness and D. Jespersen
   !       "Moderate Degree Symmetric Quadrature Rules for the Triangle"
   !       J. Inst. Math. Appl., 1975, 15(1), pp. 19-32
   !
   !       "J. Inst. Math. Appl." is now Continued as "IMA J. Appl. Math."
   !       J. Inst. Math. Appl. = Journal of the Institute of Mathematics and its Applications
   !       IMA J. Appl. Math.   = IMA Journal of Applied Mathematics
   !
   !***********************************************************************
   !
   implicit none
   integer(8) nquad, nquad_max
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)

   !     Local variables
   integer(8) i
   double precision wgt_bar, coor_bar
   double precision wgt_1, coor_1_grp_1, coor_2_grp_1
   double precision wgt_2, coor_1_grp_2, coor_2_grp_2
   double precision wgt_3, coor_1_grp_3, coor_2_grp_3
   double precision wgt_4, coor_1_grp_4, coor_2_grp_4, coor_3_grp_4


   !     Check the number of quadrature points

   nquad = 16
   if (nquad .gt. nquad_max) then
      write(*,*)
      write(*,*) "   ???"
      write(*,*) "quad_triangle: nquad > nquad_max : ",&
      &nquad, nquad_max
      write(*,*) "quad_triangle: Aborting..."
      stop
   endif



   !_____________________________________________________________________

   wgt_bar = 1.443156076777862d-1 / 2.d0
   wgt_1   = 2.852749028018549d-1 / 6.d0
   wgt_2   = 9.737549286959440d-2 / 6.d0
   wgt_3   = 3.096521116041552d-1 / 6.d0
   wgt_4   = 1.633818850466092d-1 / 1.2d1

   !_____________________________________________________________________
   coor_bar = 1.d0 / 3.d0

   coor_1_grp_1 = 4.592925882927229d-1
   coor_2_grp_1 = 8.141482341455413d-2

   coor_1_grp_2 = 5.054722831703103d-2
   coor_2_grp_2 = 8.989055433659379d-1

   coor_1_grp_3 = 1.705693077517601d-1
   coor_2_grp_3 = 6.588613844964797d-1

   coor_1_grp_4 = 7.284923929554041d-1
   coor_2_grp_4 = 2.63112829634638689d-1
   coor_3_grp_4 = 8.394777409957211d-3
   !_____________________________________________________________________

   i = 1
   xq(i) = coor_bar
   yq(i) = coor_bar
   wq(i) = wgt_bar

   !_____________________________________________________________________
   i = 2
   xq(i) = coor_1_grp_1
   yq(i) = coor_1_grp_1
   wq(i) = wgt_1

   i = 3
   xq(i) = coor_1_grp_1
   yq(i) = coor_2_grp_1
   wq(i) = wgt_1

   i = 4
   xq(i) = coor_2_grp_1
   yq(i) = coor_1_grp_1
   wq(i) = wgt_1

   !_____________________________________________________________________
   i = 5
   xq(i) = coor_1_grp_2
   yq(i) = coor_1_grp_2
   wq(i) = wgt_2

   i = 6
   xq(i) = coor_1_grp_2
   yq(i) = coor_2_grp_2
   wq(i) = wgt_2

   i = 7
   xq(i) = coor_2_grp_2
   yq(i) = coor_1_grp_2
   wq(i) = wgt_2

   !_____________________________________________________________________
   i = 8
   xq(i) = coor_1_grp_3
   yq(i) = coor_1_grp_3
   wq(i) = wgt_3

   i = 9
   xq(i) = coor_1_grp_3
   yq(i) = coor_2_grp_3
   wq(i) = wgt_3

   i = 10
   xq(i) = coor_2_grp_3
   yq(i) = coor_1_grp_3
   wq(i) = wgt_3

   !_____________________________________________________________________
   i = 11
   xq(i) = coor_1_grp_4
   yq(i) = coor_2_grp_4
   wq(i) = wgt_4

   i = 12
   xq(i) = coor_2_grp_4
   yq(i) = coor_1_grp_4
   wq(i) = wgt_4

   i = 13
   xq(i) = coor_2_grp_4
   yq(i) = coor_3_grp_4
   wq(i) = wgt_4

   i = 14
   xq(i) = coor_3_grp_4
   yq(i) = coor_2_grp_4
   wq(i) = wgt_4

   i = 15
   xq(i) = coor_3_grp_4
   yq(i) = coor_1_grp_4
   wq(i) = wgt_4

   i = 16
   xq(i) = coor_1_grp_4
   yq(i) = coor_3_grp_4
   wq(i) = wgt_4

end
