C Calculate the overlap integral of two EM modes and an AC mode
C using numerical quadrature.
C
      subroutine photoelastic_int_v2 (nval_EM, nval_AC, ival1,
     *  ival2, ival3, nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, p_tensor, beta_AC, soln_EM, soln_AC, eps_lst, 
     *  debug, overlap, basis_overlap)
c
      implicit none
      integer*8 nval_EM, nval_AC, ival1, ival2, ival3
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel), debug
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
c      complex*16 x(2,npt)
      complex*16 soln_EM(3,nnodes,nval_EM,nel)
      complex*16 soln_AC(3,nnodes,nval_AC,nel)
      complex*16 overlap, beta_AC
      complex*16 p_tensor(3,3,3,3,nb_typ_el)

c     Local variables
      integer*8 nnodes0
      parameter (nnodes0 = 6)
      double precision xel(2,nnodes0)
      complex*16 basis_overlap(3*nnodes0,3*nnodes0,3,3*nnodes0)
      complex*16 E1star, E2, Ustar, eps
      integer*8 i, j, k, l, j1, typ_e
      integer*8 iel, ind_ip, i_eq
      integer*8 jtest, ind_jp, j_eq, k_eq
      integer*8 ltest, ind_lp, l_eq
      integer*8 itrial, ui
      complex*16 eps_lst(nb_typ_el)
      complex*16 z_tmp1, ii
      double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
      double precision det_b, eps_0
c
c     NQUAD: The number of quadrature points used in each element.
      integer*8 nquad, nquad_max
      parameter (nquad_max = 16) ! Limit to P2 polynomials
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
cc      integer*8 info_curved, n_curved
      double precision ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
      complex*16 coeff

      double precision p2_p2_p2(6,6,6)
      double precision p2_p2_p2x(6,6,6), p2_p2_p2y(6,6,6)
C
C
Cf2py intent(in) nval_EM, nval_AC, ival1, ival2, ival3, nb_typ_el
Cf2py intent(in) nel, npt, nnodes, table_nod, p_tensor, beta_AC , debug
Cf2py intent(in) type_el, x, soln_EM, soln_AC, eps_lst
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_EM) nnodes, nval_EM, nel
Cf2py depend(soln_AC) nnodes, nval_AC, nel
Cf2py depend(p_tensor) nb_typ_el
C
Cf2py intent(out) overlap, basis_overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      eps_0 = 1.0d0!8.854187817d-12
      ii = cmplx(0.0d0, 1.0d0)
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "photoelastic_int_v2: problem nnodes = ", nnodes
        write(ui,*) "photoelastic_int_v2: nnodes should be equal to 6 !"
        write(ui,*) "photoelastic_int_v2: Aborting..."
        stop
      endif
C
      overlap = 0.0d0
      call quad_triangle (nquad, nquad_max, wq, xq, yq)
      if (debug .eq. 1) then
        write(ui,*) "photoelastic_int_v2: nquad, nquad_max = ",
     *              nquad, nquad_max
      endif
cccccccccccc
C Loop over elements - start
cccccccccccc
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
cccccccccc
cccccccc
c       The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
c       maps the current triangle to the reference triangle.
        do i=1,2
          do j=1,2
            mat_B(j,i) = xel(j,i+1) - xel(j,1)
          enddo
        enddo
        det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
        if (abs(det_b) .le. 1.0d-22) then  ! TEMPORARY CHANGE
cc        if (abs(det_b) .le. 1.0d-8) then
          write(*,*) '?? PE_int_v2: Determinant = 0 :', det_b
          write(*,*) "xel = ", xel
          write(*,*) 'Aborting...'
          stop
        endif
c       We also need, is the matrix mat_T of the reverse transmation 
c                (from reference to current triangle):
c       mat_T = inverse matrix of de mat_B
        mat_T(1,1) =  mat_B(2,2) / det_b
        mat_T(2,2) =  mat_B(1,1) / det_b
        mat_T(1,2) = -mat_B(1,2) / det_b
        mat_T(2,1) = -mat_B(2,1) / det_b
c       Note that if grad_i_0 is the gradient on the reference triangle, 
c       then the gradient on the actual triangle is:
c       grad_i  = Transpose(mat_T)*grad_i0
c
c       mat_T_tr = Transpose(mat_T)
        mat_T_tr(1,1) = mat_T(1,1)
        mat_T_tr(2,2) = mat_T(2,2)
        mat_T_tr(1,2) = mat_T(2,1)
        mat_T_tr(2,1) = mat_T(1,2)

      call mat_p2_p2_p2 (p2_p2_p2, det_b)
      call mat_p2_p2_p2x (p2_p2_p2x, mat_T_tr, det_b)
      call mat_p2_p2_p2y (p2_p2_p2y, mat_T_tr, det_b)

cccccccccc
        do i=1,3*nnodes
          do j=1,3*nnodes
            do k=1,3
              do l=1,3*nnodes
                basis_overlap(i,j,k,l) = 0.0d0
              enddo
            enddo
          enddo
        enddo
cccccccccc
          do itrial=1,nnodes0
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
              do jtest=1,nnodes0
                do j_eq=1,3
                  ind_jp = j_eq + 3*(jtest-1)
C                 Gradient of transverse components of basis function
                  do k_eq=1,3
                    do ltest=1,nnodes0
                      do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
                        if ( k_eq .eq. 1) then
                          z_tmp1 = p2_p2_p2x(itrial,jtest,ltest)
                        elseif ( k_eq .eq. 2) then
                          z_tmp1 = p2_p2_p2y(itrial,jtest,ltest)
                        elseif ( k_eq .eq. 3) then
                          z_tmp1 = p2_p2_p2(itrial,jtest,ltest)
                          z_tmp1 = z_tmp1 * ii * beta_AC
                        else
                          write(*,*) "--- photoelastic_int_v2: "
                          write(*,*) "k_eq has illegal value:"
                          write(*,*) "k_eq = ", k_eq
                          write(*,*) "Aborting..."
                          stop
                        endif
                        coeff = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                        eps = eps_lst(typ_e)
                        z_tmp1 = coeff * eps**2 * z_tmp1
                        basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) =
     *                          z_tmp1
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo


cc      stop


cccccccc
cccccccccc
cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do itrial=1,nnodes0
          do i_eq=1,3
            ind_ip = i_eq + 3*(itrial-1)
            E1star = conjg(soln_EM(i_eq,itrial,ival1,iel))
            do jtest=1,nnodes0
              do j_eq=1,3
                ind_jp = j_eq + 3*(jtest-1)
                E2 = soln_EM(j_eq,jtest,ival2,iel)
                do ltest=1,nnodes0
                  do l_eq=1,3
                    ind_lp = l_eq + 3*(ltest-1)
                    Ustar = conjg(soln_AC(l_eq,ltest,ival3,iel))
                    do k_eq=1,3
                      z_tmp1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                      z_tmp1 = E1star * E2 * Ustar * z_tmp1
                      overlap = overlap + z_tmp1
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
cccccccccccc
C Loop over elements - end
cccccccccccc
      enddo
C Apply scaling that sits outside of integration.
      overlap = overlap * eps_0
      if(debug .eq. 1) then
        write(*,*) "PE_int_v2: overlap"
        write(*,*) overlap
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       open (unit=26,file="Output/overlap_v2.txt")
C         write(26,*) "overlap, eps_0 = "
C         write(26,*) overlap, eps_0
C       close (unit=26)
C         open(4,file="Output/basis_overlap_v2.txt",status='unknown')
C           iel = nel
C           do itrial=1,nnodes0
C             do i_eq=1,3
C               ind_ip = i_eq + 3*(itrial-1)
C               do j_eq=1,3
C                 do k_eq=1,3
C                   do ltest=1,nnodes0
C                     do l_eq=1,3
C                       ind_lp = l_eq + 3*(ltest-1)
C                       z_tmp1  = basis_overlap(ind_ip,j_eq,k_eq,ind_lp)
C                       if (z_tmp1 .ne. 0) then
C                         write(4,*) ind_ip,j_eq,k_eq,ind_lp,
C      *                  abs(z_tmp1), z_tmp1
C                       endif
C                     enddo
C                   enddo
C                 enddo
C               enddo
C             enddo
C           enddo
C         close(4)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      end subroutine photoelastic_int_v2