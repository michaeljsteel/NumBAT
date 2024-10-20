#include "numbat_decl.h"


subroutine photoelastic_int_curvilinear_elts (nval_em_p, nval_em_s, nval_ac_u, &
   ival_p, ival_s, ival_ac, &
   n_msh_el, n_msh_pts, elnd_to_mesh, v_nd_xy, &
   n_elt_mats, el_material, p_tensor, beta_ac, soln_em_p, soln_em_s, soln_ac_u, &
   v_eps_rel, Q_PE, errco, emsg)

   use numbatmod
   use alloc

   use class_QuadIntegrator

   integer(8) n_msh_el, n_msh_pts,  n_elt_mats
   integer(8) nval_em_p, nval_em_s, nval_ac_u, ival_p, ival_s, ival_ac

   integer(8) el_material(n_msh_el)
   integer(8) elnd_to_mesh(P2_NODES_PER_EL, n_msh_el)
   double precision v_nd_xy(2, n_msh_pts)

   complex(8) soln_em_p(3, P2_NODES_PER_EL, nval_em_p, n_msh_el)
   complex(8) soln_em_s(3, P2_NODES_PER_EL, nval_em_s, n_msh_el)
   complex(8) soln_ac_u(3, P2_NODES_PER_EL, nval_ac_u, n_msh_el)

   complex(8) p_tensor(3,3,3,3,n_elt_mats)
   complex(8) beta_ac
   complex(8) v_eps_rel(n_elt_mats)
   complex(8), intent(out) :: Q_PE(nval_em_s, nval_em_p, nval_ac_u)

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !     Local variables

   double precision nds_xy(2,P2_NODES_PER_EL)

   ! Contains phi2_i phi2_j \partial_k phi2_l for each basis function
   ! Slots i,j, ,l are indexed phi_1x, phi_1y, phi_1z, phi_2x, phi_2y, ...
   !complex(8), dimension(:,:,:,:), allocatable :: basovrlp_old


   complex(8) E_s_i_star, E_p_j, Ustar_l, eps
   integer(8)  j,    typ_e
   integer(8) i_el, ind_ip, xyz_i
   integer(8) nd_j, ind_jp, xyz_j, xyz_k
   integer(8) nd_l, ind_lp, xyz_l
   integer(8) nd_i, ui
   complex(8) zt1
!   double precision mat_B(2,2), mat_T(2,2)


   !     NQUAD: The number of quadrature points used in each element.
   ! integer(8) nquad, iq

   ! Limit to P2 polynomials, in assembly this is 25?


   ! integer(8), parameter :: NQUAD_MAX = 16
   ! double precision wt_quad(NQUAD_MAX)
   ! double precision x_quad(NQUAD_MAX), y_quad(NQUAD_MAX)
   ! double precision xy_ref(2), xy_act(2), ww, det


   integer(8) n_curved
   logical is_curved

   !complex(8) coeff_1, coeff_2
   ! double precision phi_P2_ref(P2_NODES_PER_EL)
   ! double precision gradt_P2_ref(2, P2_NODES_PER_EL), gradt_P2_act(2, P2_NODES_PER_EL)


   integer(8) v_ival_p(nval_em_p), v_ival_s(nval_EM_s), v_ival_u(nval_ac_u)
   integer(8) ivs, ivp, ivu, t_ival_s, t_ival_p, t_ival_u

   type(QuadIntegrator) quadint

   complex(8), dimension(:,:,:,:), allocatable :: basis_overlap
   !double precision reldiff, bo1, bo2

   !fo2py intent(in) nval_em_p, nval_em_s, nval_ac_u
   !fo2py intent(in) ival_p, ival_s, ival_ac, n_elt_mats
   !fo2py intent(in) n_msh_el, n_msh_pts, P2_NODES_PER_EL, elnd_to_mesh, p_tensor, beta_ac , debug
   !fo2py intent(in) el_material, x, soln_em_p, soln_em_s, soln_ac_u, v_eps_rel

   !f2py depend(elnd_to_mesh) P2_NODES_PER_EL, n_msh_el
   !f2py depend(el_material) n_msh_pts
   !f2py depend(x) n_msh_pts
   !f2py depend(soln_em_p) P2_NODES_PER_EL, nval_em_p, n_msh_el
   !f2py depend(soln_em_s) P2_NODES_PER_EL, nval_em_s, n_msh_el
   !f2py depend(soln_ac_u) P2_NODES_PER_EL, nval_ac_u, n_msh_el
   !f2py depend(p_tensor) n_elt_mats
   !f2py depend(v_eps_rel) n_elt_mats

   !fo2py intent(out) Q_PE
   !fo2py intent(out) errco
   !fo2py intent(out) emsg



   ui = stdout


   Q_PE = D_ZERO

   ! build arrays holding which modes to be calculated
   call fill_ival_arrays(v_ival_p, nval_em_p, ival_p)
   call fill_ival_arrays(v_ival_s, nval_EM_s, ival_s)
   call fill_ival_arrays(v_ival_u, nval_ac_u, ival_ac)


   call quadint%setup_reference_quadratures()

   !call complex_alloc_4d(basovrlp_old, 3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, &
   !   3_8, 3*P2_NODES_PER_EL, 'basovrlp_old', errco, emsg)

   call complex_alloc_4d(basis_overlap, 3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, &
      3_8, 3*P2_NODES_PER_EL, 'basovrlp_old', errco, emsg)

   Q_PE = D_ZERO


   do i_el=1,n_msh_el
      typ_e = el_material(i_el)
      eps = v_eps_rel(typ_e)

      ! find positions of all the P2 nodes for this elt
      do j=1,P2_NODES_PER_EL
         nds_xy(:, j) = v_nd_xy(:,  elnd_to_mesh(j,i_el))
      enddo

      is_curved = log_is_curved_elem_tri (P2_NODES_PER_EL, nds_xy)
      if (is_curved) then
         n_curved = n_curved + 1
      endif


      call make_P2_overlaps_i_j_dk_l(i_el, beta_ac, typ_e, is_curved, &
         eps, p_tensor, n_elt_mats, nds_xy, quadint, basis_overlap, errco, emsg)
      RETONERROR(errco)


      ! basovrlp_old = D_ZERO

      ! ! For each quadrature point evaluate Q_PE of Lagrange polynomials
      ! ! or derivative of Lagrange polynomials
      ! do iq=1,nquad
      !    xy_ref(1) = x_quad(iq)
      !    xy_ref(2) = y_quad(iq)
      !    ww = wt_quad(iq)
      !    !  xy_ref   = coordinate on the reference triangle
      !    !  xy_act = coordinate on the actual triangle

      !    !  phi_P2_ref = values of Lagrange polynomials (1-6) at each local node.
      !    !  gradt_P2_ref = gradient on the reference triangle (P2 element)
      !    call phi2_2d_mat(xy_ref, phi_P2_ref, gradt_P2_ref)

      !    if (.not. is_curved) then ! Rectilinear element
      !       call jacobian_p1_2d(xy_ref, nds_xy, P2_NODES_PER_EL, xy_act, det, mat_B, mat_T, errco, emsg)
      !    else ! Isoparametric element !  fixed 2024/6/12
      !       call jacobian_p2_2d(nds_xy, P2_NODES_PER_EL, phi_P2_ref, gradt_P2_ref, xy_act, det, mat_B, mat_T, errco, emsg)
      !    endif

      !    if(abs(det) .lt. 1.0d-20) then
      !       write(*,*)
      !       write(*,*) "   ???"
      !       write(*,*) "PE_int: det = 0 : i_el, det = ", i_el, det
      !       write(*,*) "PE_int: Aborting..."
      !       stop
      !    endif

      !    ! grad_i  = gradient on the actual triangle
      !    ! grad_i  = Transpose(mat_T)*grad_i0
      !    ! Calculation of the matrix-matrix product:
      !    call DGEMM('Transpose','N', 2, 6, 2, D_ONE, mat_T, 2, gradt_P2_ref, 2, D_ZERO, gradt_P2_act, 2)
      !    coeff_1 = ww * abs(det)



      !    ! Calculate Q_PE of basis functions at quadrature point,
      !    ! which is a superposition of P2 polynomials for each function (fi_eld).
      !    do nd_i=1,P2_NODES_PER_EL
      !       do xyz_i=1,3                        ! x/y/z component
      !          ind_ip = xyz_i + 3*(nd_i-1)

      !          do nd_j=1,P2_NODES_PER_EL
      !             do xyz_j=1,3                      ! x/y/z component
      !                ind_jp = xyz_j + 3*(nd_j-1)

      !                ! Gradient of transverse components of basis function
      !                do xyz_k=1,2                    ! x/y component
      !                   do nd_l=1,P2_NODES_PER_EL
      !                      do xyz_l=1,3
      !                         ind_lp = xyz_l + 3*(nd_l-1)
      !                         zt1 = phi_P2_ref(nd_i) * phi_P2_ref(nd_j) * gradt_P2_act(xyz_k,nd_l)
      !                         coeff_2 = p_tensor(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)

      !                         zt1 = coeff_1 * coeff_2 * eps**2 * zt1


      !                         basovrlp_old(ind_ip,ind_jp,xyz_k,ind_lp) = basovrlp_old(ind_ip,ind_jp,xyz_k,ind_lp) + zt1


      !                         ! if (i_el .eq. 1 .and.  &
      !                         !    ind_ip .eq. 1 .and. ind_jp .eq. 1 .and. &
      !                         !    xyz_k .eq. 1 .and.  ind_lp .eq. 2 ) then
      !                         !    write(*,'("oldinner ", i3, " (", e12.4,e12.4, ") ", "(", e12.4,e12.4, ")")') &
      !                         !       iq, zt1, basovrlp_old(ind_ip, ind_jp, xyz_k, ind_lp)
      !                         ! endif



      !                      enddo
      !                   enddo
      !                enddo

      !                ! Gradient of longitudinal components of basis function,
      !                !  which is i*beta*phi because field is assumed to be of form e^{i*beta*z} phi.
      !                xyz_k=3                        ! z component
      !                do nd_l=1,P2_NODES_PER_EL
      !                   do xyz_l=1,3
      !                      ind_lp = xyz_l + 3*(nd_l-1)
      !                      zt1 = phi_P2_ref(nd_i) * phi_P2_ref(nd_j) * phi_P2_ref(nd_l) * (-C_IM_ONE* beta_ac)
      !                      coeff_2 = p_tensor(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)
      !                      zt1 = coeff_1 * coeff_2 * eps**2 * zt1


      !                      basovrlp_old(ind_ip,ind_jp,xyz_k,ind_lp) = basovrlp_old(ind_ip,ind_jp,xyz_k,ind_lp) + zt1

      !                      ! if (i_el .eq. 1 .and. &
      !                      !    ind_ip .eq. 1 .and. ind_jp .eq. 1 .and. &
      !                      !    xyz_k .eq. 1 .and.  ind_lp .eq. 2 ) then
      !                      !    write(*,'("oldinner ", i3, " (", e12.4,e12.4, ") ", "(", e12.4,e12.4, ")")') &
      !                      !       iq, zt1, basovrlp_old(ind_ip, ind_jp, xyz_k, ind_lp)
      !                      ! endif



      !                   enddo
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    enddo
      ! enddo

      ! ! if(i_el .eq. 1) then
      ! !    do xyz_i=1,3*P2_NODES_PER_EL
      ! !       do xyz_j=1,3*P2_NODES_PER_EL
      ! !          do xyz_k=1,3
      ! !             do xyz_l=1,3*P2_NODES_PER_EL
      ! !                bo1 = basovrlp_old(xyz_i, xyz_j, xyz_k, xyz_l)
      ! !                bo2 = basis_overlap(xyz_i, xyz_j, xyz_k, xyz_l)
      ! !                !reldiff = dabs( (bo1-bo2)/bo1)
      ! !                reldiff = dabs(bo1-bo2)
      ! !                if (reldiff .gt. 1e-20) then
      ! !                   write(*, '(A,i4,i4,i4,i4,e13.3,e13.3,e13.3,e13.3,e13.3)') 'bdiff ', &
      ! !                      xyz_i, xyz_j, xyz_k, xyz_l, bo1, bo2, reldiff
      ! !                endif
      ! !             enddo
      ! !          enddo
      ! !       enddo
      ! !    enddo

      ! ! endif


      ! Having calculated Q_PE of basis functions on element
      ! now multiply by specific field values for modes of interest.
      !


      do nd_i=1,P2_NODES_PER_EL          ! nodes and components of the Stokes field
         do xyz_i=1,3
            ind_ip = xyz_i + 3*(nd_i-1)

            do ivs = 1,nval_em_s          ! for the requested Stokes modes
               t_ival_s = v_ival_s(ivs)
               if (t_ival_s .eq. 0) exit

               ! (E_S^*)_i
               E_s_i_star = conjg(soln_em_s(xyz_i, nd_i, t_ival_s, i_el))

               do nd_j=1,P2_NODES_PER_EL  ! nodes and components of the pump field
                  do xyz_j=1,3
                     ind_jp = xyz_j + 3*(nd_j-1)

                     do ivp = 1,nval_em_p   ! for the requested pump modes
                        t_ival_p = v_ival_p(ivp)
                        if (t_ival_p .eq. 0) exit

                        ! (E_P^)_j
                        E_p_j = soln_em_p(xyz_j, nd_j, t_ival_p, i_el)


                        do nd_l=1,P2_NODES_PER_EL  ! nodes and components of the AC field
                           do xyz_l=1,3
                              ind_lp = xyz_l + 3*(nd_l-1)

                              do ivu = 1,nval_ac_u
                                 t_ival_u = v_ival_u(ivu)
                                 if (t_ival_u .eq. 0) exit

                                 Ustar_l = conjg(soln_ac_u(xyz_l,nd_l,t_ival_u,i_el))
                                 zt1 = E_s_i_star * E_p_j * Ustar_l

                                 do xyz_k=1,3
                                    Q_PE(t_ival_s,t_ival_p,t_ival_u) = Q_PE(t_ival_s,t_ival_p,t_ival_u) &
                                       + zt1 * basis_overlap(ind_ip,ind_jp,xyz_k,ind_lp)
                                 enddo

                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   Q_PE = - Q_PE * SI_EPS_0


end subroutine photoelastic_int_curvilinear_elts

