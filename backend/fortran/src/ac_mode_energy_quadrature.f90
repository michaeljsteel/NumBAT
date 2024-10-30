! Calculate the elastic energy v_energy integral of an AC mode with itself using
! numerical quadrature.

! E = 2 \Omega^2 \int dxdy \rho |u|^2

subroutine AC_mode_energy_quadrature (n_modes, n_msh_el, n_msh_pts, v_nd_xy, &
   elnd_to_mshpt, n_elt_mats, v_el_material,  &
   rho, Omega_AC, soln_ac_u, debug, v_energy)

   use numbatmod
   integer(8) n_modes, ival
   integer(8) n_msh_el, n_msh_pts, n_elt_mats
   integer(8) v_el_material(n_msh_el), debug
   integer(8) elnd_to_mshpt(P2_NODES_PER_EL,n_msh_el)
   double precision v_nd_xy(2,n_msh_pts)
   complex(8) soln_ac_u(3,P2_NODES_PER_EL,n_modes,n_msh_el)
   complex(8) Omega_AC(n_modes)
   complex(8), dimension(n_modes) :: v_energy
   complex(8) rho(n_elt_mats)

   integer(8) errco
   character(len=EMSG_LENGTH) emsg

!     Local variables

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) bas_ovrlp(P2_NODES_PER_EL,P2_NODES_PER_EL)
   complex(8) U, Ustar
   integer(8) i, j, j1, typ_e
   integer(8) i_el, k_eq
   integer(8) itrial, jtest, ui
   double precision mat_B(2,2), mat_T(2,2)
!
!     NQUAD: The number of quadrature points used in each element.
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   logical  is_curved
   integer(8) n_curved
   complex(8) coeff_1
   double precision phi2_list(6), grad2_mat0(2,6)
!
!
!f2py intent(in) n_modes, n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mshpt
!f2py intent(in) v_el_material, x, n_elt_mats, rho
!f2py intent(in) soln_ac_u, debug, Omega_AC
!
!f2py depend(elnd_to_mshpt) P2_NODES_PER_EL, n_msh_el
!f2py depend(v_el_material) n_msh_pts
!f2py depend(x) n_msh_pts
!f2py depend(soln_ac_u) P2_NODES_PER_EL, n_modes, n_msh_el
!f2py depend(Omega_AC) n_modes
!f2py depend(rho) n_elt_mats
!
!f2py intent(out) v_energy
!

   ui = stdout

   call quad_triangle (nquad, nquad_max, wq, xq, yq)
   if (debug .eq. 1) then
      write(ui,*) "AC_mode_elastic_energy_int: nquad, nquad_max = ",&
      &nquad, nquad_max
   endif

   v_energy = D_ZERO


   do i_el=1,n_msh_el
      typ_e = v_el_material(i_el)
      do j=1,P2_NODES_PER_EL
         j1 = elnd_to_mshpt(j,i_el)
         nod_el_p(j) = j1
         xel(:,j) = v_nd_xy(:,j1)
      enddo

      is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, xel)
      if (is_curved ) then
         n_curved = n_curved + 1
      endif


      do j=1,P2_NODES_PER_EL
         do i=1,P2_NODES_PER_EL
            bas_ovrlp(i,j) = 0.0d0
         enddo
      enddo

! For each quadrature point evaluate v_energy of Lagrange polynomials
! or derivative of Lagrange polynomials
      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
!         xx   = coordinate on the reference triangle
!         xx_g = coordinate on the actual triangle
!         phi2_list = values of Lagrange polynomials (1-6) at each local node.
!         grad2_mat0 = gradient on the reference triangle (P2 element)
         call phi2_2d_mat(xx, phi2_list, grad2_mat0)
!
         if (.not. is_curved ) then
!           Rectilinear element
            call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL,&
            &xx_g, det, mat_B, mat_T, errco, emsg)
         else
!           Isoparametric element, 2024-06-13 fixed version
            call jacobian_p2_2d(xel, P2_NODES_PER_EL, phi2_list,&
            &grad2_mat0, xx_g, det, mat_B, mat_T, errco, emsg)
         endif
         if(abs(det) .lt. 1.0d-20) then
            write(*,*)
            write(*,*) "   ???"
            write(*,*) "AC_m_el_en_int: det = 0 : i_el, det =", i_el, det
            write(*,*) "AC_m_el_en_int: Aborting..."
            stop
         endif

         coeff_1 = ww * abs(det)

! Calculate v_energy of basis functions at quadrature point,
! which is a superposition of P2 polynomials for each function (field).
         do itrial=1,P2_NODES_PER_EL
            do jtest=1,P2_NODES_PER_EL
               bas_ovrlp(itrial,jtest) =&
               &bas_ovrlp(itrial,jtest) +&
               &coeff_1 * phi2_list(itrial) * phi2_list(jtest)
            enddo
         enddo
      enddo

! Having calculated v_energy of basis functions on element
! now multiply by specific field values for modes of interest.
      do ival=1,n_modes
         do itrial=1,P2_NODES_PER_EL
            do jtest=1,P2_NODES_PER_EL
               do k_eq=1,3
                  Ustar = conjg(soln_ac_u(k_eq,itrial,ival,i_el))
                  U = soln_ac_u(k_eq,jtest,ival,i_el)
                  v_energy(ival) = v_energy(ival) +&
                  &rho(typ_e) * Ustar * U * bas_ovrlp(itrial,jtest)
               enddo
            enddo
         enddo
      enddo

   enddo

   v_energy = 2.0 * Omega_AC**2 * v_energy

end subroutine AC_mode_energy_quadrature
