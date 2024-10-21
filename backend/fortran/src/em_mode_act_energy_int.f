! Calculate the energy (not power) overlap integral of an EM mode with itself using
! numerical quadrature.
!
subroutine em_mode_act_energy_int (nval, nel, npt,&
&elnd_to_mesh, type_el, nb_typ_el, n_lst,&
&x, soln_EM, overlap)
!
   use numbatmod
   integer(8) nval, nel, npt
   integer(8) elnd_to_mesh(P2_NODES_PER_EL,nel), nb_typ_el
   integer(8) type_el(nel)
   double precision x(2,npt)
   complex(8) soln_EM(3,P2_NODES_PER_EL+7,nval,nel)
   complex(8) n_lst(nb_typ_el), eps_lst(nb_typ_el)
   complex(8), dimension(nval) :: overlap
   integer(8) errco
   character(len=EMSG_LENGTH) emsg

!     Local variables
   integer(8) nod_el_p(P2_NODES_PER_EL)
   complex(8) basis_overlap(P2_NODES_PER_EL)
   integer(8) i, j, j1, typ_e
   integer(8) iel, ival
   integer(8) itrial, i_eq
   logical  is_curved
   integer(8) n_curved, debug, ui
   double precision xel(2,P2_NODES_PER_EL)
   double precision phi2_list(6), grad2_mat0(2,6)
   double precision grad2_mat(2,6)

   complex(8) coeff_1
   complex(8) E, Estar
!
!     NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   double precision mat_B(2,2), mat_T(2,2)
!
!
!f2py intent(in) nval, nel, npt
!f2py intent(in) P2_NODES_PER_EL, elnd_to_mesh, type_el, nb_typ_el
!f2py intent(in) x, soln_EM, n_lst
!
!f2py depend(elnd_to_mesh) P2_NODES_PER_EL, nel
!f2py depend(x) npt
!f2py depend(soln_EM) P2_NODES_PER_EL, nval, nel
!f2py depend(n_lst) nb_typ_el
!f2py depend(type_el) nel
!
!f2py intent(out) overlap
!
!
   !!!!!!!!!!!!!!!!!!!!!!!!  Start Program  !!!!!!!!!!!!!!!!!!!!!!!!
!
   ui = stdout
   debug = 0

!
   if ( P2_NODES_PER_EL .ne. 6 ) then
      write(ui,*) "EM_mode_E_energy_int: problem P2_NODES_PER_EL = ", P2_NODES_PER_EL
      write(ui,*)"EM_mode_E_energy_int: P2_NODES_PER_EL should be equal to 14!"
      write(ui,*) "EM_mode_E_energy_int: Aborting..."
      stop
   endif
!
   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "EM_mode_E_energy_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

!     Calculate permittivity
   do i = 1, int(nb_typ_el)
      eps_lst(i) = n_lst(i)**2
   enddo
!
   do i=1,nval
      overlap(i) = 0.0d0
   enddo
!
!      n_curved = 0
   do iel=1,nel
      typ_e = type_el(iel)
      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mesh(j,iel)
         nod_el_p(j) = j1
         xel(1,j) = x(1,j1)
         xel(2,j) = x(2,j1)
      enddo
      is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, xel)
      if (is_curved) then
         n_curved = n_curved + 1
      endif
!ccccccccc
      do i=1,P2_NODES_PER_EL
         basis_overlap(i) = 0.0d0
      enddo
!ccccccccc
      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
!         xx   = coordinate on the reference triangle
!         xx_g = coordinate on the actual triangle
!         We will also need the gradients of the P1 element
!          grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)
!
         if (.not. is_curved) then
!           Rectilinear element
            call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL,&
            &xx_g, det, mat_B, mat_T, errco, emsg)
!            if (det .le. 0) then
            if (det .le. 0 .and. debug .eq. 2 .and. iq .eq. 1) then
               write(*,*) "   !!!"
               write(*,*) "EM_m_e_en_int: det <= 0: iel, det ", iel, det
            endif
         else
!           Isoparametric element , 2024-06-13 fixed version
            call jacobian_p2_2d(xel, P2_NODES_PER_EL, phi2_list,&
            &grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif
         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "EM_m_e_en_int: det = 0 : iel, det = ", iel, det
            write(*,*) "EM_m_e_en_int: Aborting..."
            stop
         endif
!          grad_i  = gradient on the actual triangle
!          grad_i  = Transpose(mat_T)*grad_i0
!          Calculation of the matrix-matrix product:
         call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2,&
         &grad2_mat0, 2, D_ZERO, grad2_mat, 2)
         coeff_1 = ww * abs(det)
! Calculate overlap of basis functions at quadrature point,
! which is a superposition of P2 polynomials for each function (field).
         do itrial=1,P2_NODES_PER_EL
            basis_overlap(itrial) = basis_overlap(itrial) +&
            &coeff_1 * phi2_list(itrial) * phi2_list(itrial)
         enddo
      enddo
!ccccccccc
! Having calculated overlap of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,nval
         do itrial=1,P2_NODES_PER_EL
            do i_eq=1,3
               Estar = conjg(soln_EM(i_eq,itrial,ival,iel))
               E = soln_EM(i_eq,itrial,ival,iel)
               overlap(ival) = overlap(ival) +&
               &eps_lst(typ_e) * Estar * E * basis_overlap(itrial)
            enddo
         enddo
      enddo
!ccccccccccc
! Loop over elements - end
!ccccccccccc
   enddo
! Multiply through prefactor
   do i=1,nval
      overlap(i) = 2.0 * overlap(i) * SI_EPS_0
   enddo
!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
end subroutine EM_mode_E_energy_int
