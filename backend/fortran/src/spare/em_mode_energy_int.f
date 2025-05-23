C Calculate the overlap integral of an EM mode with itself using
C numerical quadrature.
C
      subroutine EM_mode_energy_int (k_0, nval, nel, npt,
     *  nnodes, elnd_to_mesh,
     *  x, betas, soln_k1, overlap)
c
C     k_0 = 2 pi / lambda, where lambda in meters.
C
      use numbatmod
      integer(8) nval, nel, npt, nnodes
      integer(8) elnd_to_mesh(nnodes,nel)
      double precision x(2,npt)
      complex(8) soln_k1(3,nnodes+7,nval,nel)
      complex(8) beta1
      complex(8) betas(nval)
      complex(8), dimension(nval) :: overlap
      double precision k_0

c     Local variables

      integer(8) nod_el_p(nnodes_0)
      complex(8) sol_el_1(2*nnodes_0+10), sol_el_2(2*nnodes_0)
      complex(8) vec_1(2*nnodes_0)
      complex(8) basis_overlap(2*nnodes_0,2*nnodes_0+10)
      integer(8) i, j, j1
      integer(8) iel, ival
      integer(8) jtest, ind_jp, j_eq
      integer(8) itrial, ind_ip, i_eq
      integer(8) n_curved, debug, ui
      logical is_curved

      double precision xel(2,nnodes_0)
      double precision phi2_list(6), grad2_mat0(2,6)
      double precision grad2_mat(2,6)
      double precision phi3_list(10), grad3_mat0(2,10)
      double precision grad3_mat(2,10)
      double precision vec_phi_j(2), vec_phi_i(2)

      complex(8) z_tmp1, z_tmp2, coeff_1
c
c     NQUAD: The number of quadrature points used in each element.
      integer(8) nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      double precision mat_B(2,2), mat_T(2,2)
C
C
Cf2py intent(in) k_0, nval, nel, npt
Cf2py intent(in) nnodes, elnd_to_mesh
Cf2py intent(in) x, betas, soln_k1
C
Cf2py depend(elnd_to_mesh) nnodes, nel
Cf2py depend(x) npt
Cf2py depend(betas) nval
Cf2py depend(soln_k1) nnodes, nval, nel
C
Cf2py intent(out) overlap
C
C
 !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
C
      ui = stdout

      debug = 0
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "EM_mode_energy_int: problem nnodes = ", nnodes
        write(ui,*) "EM_mode_energy_int: nnodes should be equal to 14 !"
        write(ui,*) "EM_mode_energy_int: Aborting..."
        stop
      endif
c
      call quad_triangle (nquad, nquad_max, wq, xq, yq)
      if (debug .eq. 1) then
        write(ui,*) "EM_mode_energy_int: nquad, nquad_max = ",
     *              nquad, nquad_max
      endif
C
      do i=1,nval
        overlap(i) = 0.0d0
      enddo
c
c      n_curved = 0
      do iel=1,nel
        do j=1,nnodes
          j1 = elnd_to_mesh(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        is_curved = log_is_curved_elem_tri (nnodes, xel)
        if (is_curved) then
          n_curved = n_curved + 1
        endif
cccccccccc
        do i=1,2*nnodes
          do j=1,2*nnodes+10
            basis_overlap(i,j) = 0.0d0
          enddo
        enddo
cccccccccc
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
c         xx   = coordinate on the reference triangle
c         xx_g = coordinate on the actual triangle
c         We will also need the gradients of the P1 element
c          grad2_mat0 = gradient on the reference triangle (P2 element)
           call phi2_2d_mat(xx, phi2_list, grad2_mat0)
c          grad3_mat0 = gradient on the reference triangle (P3 element)
           call phi3_2d_mat(xx, phi3_list, grad3_mat0)
c
          if (.not. is_curved ) then
c           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,
     *               xx_g, det, mat_B, mat_T)
c            if (det .le. 0) then
            if (det .le. 0 .and. debug .eq. 2 .and. iq .eq. 1) then
              write(*,*) "   !!!"
              write(*,*) "EM_m_en_int: det <= 0: iel, det ", iel, det
            endif
          else
c           Isoparametric element, 2024-06-13 fixed version
            call jacobian_p2_2d(xel, nnodes, phi2_list,
     *               grad2_mat0, xx_g, det, mat_B, mat_T)
          endif
           if(abs(det) .lt. 1.0d-20) then
             write(*,*)
             write(*,*) "   ???"
             write(*,*) "EM_m_en_int: det = 0 : iel, det = ", iel, det
             write(*,*) "EM_m_en_int: Aborting..."
             stop
           endif
c          grad_i  = gradient on the actual triangle
c          grad_i  = Transpose(mat_T)*grad_i0
c          Calculation of the matrix-matrix product:
          call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,
     *      grad2_mat0, 2, D_ZERO, grad2_mat, 2)
          call DGEMM('Transpose','N', 2, 10, 2, D_ONE, mat_T, 2,
     *      grad3_mat0, 2, D_ZERO, grad3_mat, 2)
          coeff_1 = ww * abs(det)
          do itrial=1,nnodes_0
            do i_eq=1,2
              ind_ip = i_eq + 2*(itrial-1)
c             Determine the basis vector
              do i=1,2
                vec_phi_i(i) = 0.0d0
              enddo
              vec_phi_i(i_eq) = phi2_list(itrial)
              do jtest=1,nnodes_0
                do j_eq=1,2
                  ind_jp = j_eq + 2*(jtest-1)
c                 Determine the basis vector
                  do i=1,2
                    vec_phi_j(i) = 0.0d0
                  enddo
                  vec_phi_j(j_eq) = phi2_list(jtest)
                  z_tmp1 = vec_phi_i(1)*vec_phi_j(1) +
     *                        vec_phi_i(2)*vec_phi_j(2)
            z_tmp1 = coeff_1 * z_tmp1 / (k_0 * SI_C_SPEED * SI_MU_0)
                  basis_overlap(ind_ip,ind_jp) =
     *              basis_overlap(ind_ip,ind_jp) + z_tmp1
                enddo
              enddo
              do jtest=1,10
                j_eq = 3
                ind_jp = jtest + 2*nnodes_0
c               Determine the basis vector
                do i=1,2
                  vec_phi_j(i) = -grad3_mat(i,jtest)
                enddo
                z_tmp1 = vec_phi_i(1)*vec_phi_j(1) +
     *                        vec_phi_i(2)*vec_phi_j(2)
                z_tmp1 = coeff_1 * z_tmp1 / (k_0 * SI_C_SPEED * SI_MU_0)
                basis_overlap(ind_ip,ind_jp) =
     *            basis_overlap(ind_ip,ind_jp) + z_tmp1
              enddo
            enddo
          enddo
        enddo
cccccccccc
        do ival=1,nval
          beta1 = betas(ival)
          do i=1,nnodes
            do j=1,2
c             The 2 transverse components of the mode ival
              ind_ip = j + 2*(i-1)
              z_tmp1 = soln_k1(j,i,ival,iel)
c             sol_el_2 : E-field
              sol_el_2(ind_ip) = z_tmp1
            enddo
          enddo
cccccccccccc
          do i=1,nnodes
            do j=1,2
c           The 2 transverse components of the mode jval
              ind_jp = j + 2*(i-1)
              z_tmp1 = soln_k1(j,i,ival,iel)
c             sol_el_1 : H-field
              sol_el_1(ind_jp) = z_tmp1 * beta1
            enddo
          enddo
cccccccccccc
          do i=1,3
c         The longitudinal component at the vertices (P3 elements)
            ind_jp = i + 2*nnodes
            z_tmp1 = soln_k1(3,i,ival,iel)
c           sol_el_1 : H-field
            sol_el_1(ind_jp) = z_tmp1 * beta1
          enddo
          do i=nnodes+1,13
c         The longitudinal component at the edge nodes and interior node (P3 elements)
            ind_jp = i + 2*nnodes - nnodes + 3
            z_tmp1 = soln_k1(3,i,ival,iel)
            sol_el_1(ind_jp) = z_tmp1 * beta1
          enddo
cccccccccccc
c       Matrix-Vector product
          do i=1,2*nnodes
            vec_1(i) = 0.0d0
            do j=1,2*nnodes+10
              z_tmp1 = sol_el_1(j)
              z_tmp2 = basis_overlap(i,j)
              vec_1(i) = vec_1(i) + z_tmp1 * z_tmp2
            enddo
          enddo
cccccccccccc
c       Scalar product
          z_tmp1 = 0.0d0
          do i=1,2*nnodes
            z_tmp1 = vec_1(i) * sol_el_2(i)
            overlap(ival) = overlap(ival)
     *          + z_tmp1
          enddo
        enddo
      enddo
C
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
      end subroutine EM_mode_energy_int
