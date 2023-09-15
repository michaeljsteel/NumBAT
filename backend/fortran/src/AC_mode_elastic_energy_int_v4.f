C Calculate the elastic energy overlap integral of an AC mode with itself using
C analytic expressions for basis function overlaps on linear elements.  
C
      subroutine AC_mode_elastic_energy_int_v4 (nval, 
     *  nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, rho, Omega_AC, soln_AC,
     *  overlap)
c
      implicit none
      integer*8 nval, ival
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
      complex*16 soln_AC(3,nnodes,nval,nel)
      complex*16 Omega_AC(nval)
      complex*16, dimension(nval) :: overlap
      complex*16 rho(nb_typ_el)

c     Local variables
      integer*8 nnodes0
      parameter (nnodes0 = 6)
      integer*8 nod_el_p(nnodes0)
      double precision xel(2,nnodes0)
      complex*16 basis_overlap(nnodes0,nnodes0)
      complex*16 U, Ustar
      integer*8 i, j, j1, typ_e
      integer*8 iel, k_eq
      integer*8 itrial, jtest, ui
      complex*16 rho_el, z_tmp1

C
C
Cf2py intent(in) nval, nel, npt, nnodes, table_nod
Cf2py intent(in) type_el, x, nb_typ_el, rho
Cf2py intent(in) soln_AC, Omega_AC
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_AC) nnodes, nval, nel
Cf2py depend(Omega_AC) nval
Cf2py depend(rho) nb_typ_el
C
Cf2py intent(out) overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "AC_mode_elastic_energy_int_v4: problem nnodes = ", 
     *              nnodes
        write(ui,*) " --------- nnodes should be equal to 6 !"
        write(ui,*) "AC_mode_elastic_energy_int_v4: Aborting..."
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
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        rho_el = rho(typ_e)
C
        call mat_el_energy_rho (xel, rho_el, basis_overlap)

cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do ival=1,nval
          do itrial=1,nnodes0
            do jtest=1,nnodes0
              do k_eq=1,3
                Ustar = conjg(soln_AC(k_eq,itrial,ival,iel))
                U = soln_AC(k_eq,jtest,ival,iel)
                z_tmp1 = basis_overlap(itrial,jtest)
                overlap(ival)=overlap(ival)+ Ustar*U*z_tmp1
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
        overlap(i) = 2.0 * Omega_AC(i)**2 * overlap(i)
      enddo

C        open (unit=26,file="Output/overlap_elastic_v4.txt")
C        do i=1,nval
C          write(26,*) i, Omega_AC(i), abs(overlap(i)), 
C      *              overlap(i)
C        enddo
C        do i=1,nval
C          write(26,*) abs(overlap(i))
C        enddo
C        close (unit=26)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      end subroutine AC_mode_elastic_energy_int_v4
