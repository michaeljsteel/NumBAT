C Calculate the overlap integral of an AC mode with itself using
C analytic expressions for basis function overlaps on linear elements.
C
      subroutine AC_mode_power_int_v4 (nval,
     *  nel, npt, nnodes, elnd_to_mesh, type_el, x,
     *  nb_typ_el, c_tensor, beta_AC, Omega_AC, soln_AC,
     *  overlap)
c
      use numbatmod
      integer(8) nval, ival
      integer(8) nel, npt, nnodes, nb_typ_el
      integer(8) type_el(nel)
      integer(8) elnd_to_mesh(nnodes,nel)
      double precision x(2,npt)
c      complex(8) x(2,npt)
      complex(8) soln_AC(3,nnodes,nval,nel)
      complex(8) Omega_AC(nval)
      complex(8) beta_AC
      complex(8), dimension(nval) :: overlap
      complex(8) c_tensor(6,6,nb_typ_el)

c     Local variables

      integer(8) nod_el_p(P2_NODES_PER_EL)
      double precision xel(2,P2_NODES_PER_EL)
      complex(8) basis_overlap(3*P2_NODES_PER_EL,3*P2_NODES_PER_EL)
      complex(8) U, Ustar
      integer(8) i, j, j1, typ_e
      integer(8) iel, ind_ip, i_eq
      integer(8) ltest, ind_lp, l_eq
      integer(8) itrial, ui

      complex(8) z_tmp1, c_tensor_el(6,6)

C
C
Cf2py intent(in) nval, nel, npt, nnodes, elnd_to_mesh
Cf2py intent(in) type_el, x, nb_typ_el, c_tensor, beta_AC
Cf2py intent(in) soln_AC, Omega_AC
C
Cf2py depend(elnd_to_mesh) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_AC) nnodes, nval, nel
Cf2py depend(c_tensor) nb_typ_el
Cf2py depend(Omega_AC) nval
C
Cf2py intent(out) overlap
C
C
 !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
C
      ui = stdout
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "AC_mode_power_int_v4: problem nnodes = ",
     *              nnodes
        write(ui,*) " --------- nnodes should be equal to 6 !"
        write(ui,*) "AC_mode_power_int_v4: Aborting..."
        stop
      endif
C
      do i=1,nval
        overlap(i) = 0.0d0
      enddo
C
cccccccccccc
C Loop over elements - start
cccccccccccc
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = elnd_to_mesh(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        do j=1,6
          do i=1,6
            c_tensor_el(i,j) = c_tensor(i,j,typ_e)
          enddo
        enddo
        call mat_el_energy (xel, beta_AC, c_tensor_el,
     *  basis_overlap)

cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do ival=1,nval
          do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
              Ustar = conjg(soln_AC(i_eq,itrial,ival,iel))
              do ltest=1,P2_NODES_PER_EL
                do l_eq=1,3
                  ind_lp = l_eq + 3*(ltest-1)
                  U = soln_AC(l_eq,ltest,ival,iel)
                  z_tmp1 = basis_overlap(ind_ip,ind_lp)
                  overlap(ival) = overlap(ival) + Ustar * U * z_tmp1
                enddo
              enddo
            enddo
          enddo
        enddo
cccccccccccc
C Loop over elements - end
cccccccccccc
      enddo
C Multiply through prefactor
      do i=1,nval
        overlap(i) = 2.0d0 * C_IM_ONE* Omega_AC(i) * overlap(i)
      enddo

C       open (unit=26,file="Output/overlap_v4.txt")
C       do i=1,nval
C         write(26,*) i, Omega_AC(i), abs(overlap(i)),
C      *              overlap(i)
C       enddo
C       do i=1,nval
C         write(26,*) abs(overlap(i))
C       enddo
C       close (unit=26)
C
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
      end subroutine AC_mode_power_int_v4
