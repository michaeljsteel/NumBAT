#include "numbat_decl.h"

! Calculate the Q_PE integral of two EM modes and an AC mode
! using numerical quadrature over the current triangle
! defined by nds_xy.

! calculates PE style mode overlaps accounting for polarisation
! \int eps_r^2 pe(i,j,k,l) phi2_ir \phi2_js \partial_k \phi_lt  dxdy
! i,j,l run over P2 nodes
! r,s,t run over components x,y,z
!
! k runs over components x,y,z
! \partial k acts as \partial x, \partial y, or \partial z = -i beta
! see storage layout below

subroutine make_P2_overlaps_i_j_dk_l_analytic(i_el, beta_ac, typ_e,  &
   eps, p_tensor, n_elt_mats, nds_xy, bas_ovrlp)

   use numbatmod

   complex(8) beta_ac, eps
   integer(8) n_elt_mats, typ_e
   double precision nds_xy(2,P2_NODES_PER_EL)
   complex(8) p_tensor(3,3,3,3,n_elt_mats)

   ! Contains phi2_i phi2_j \partial_k phi2_l for each basis function
   ! Slots i,j, ,l are indexed phi_1x, phi_1y, phi_1z, phi_2x, phi_2y, ...
   complex(8), dimension(3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, &
      3, 3*P2_NODES_PER_EL) :: bas_ovrlp

   !locals

   integer(8) nd_i, nd_j, nd_l, xyz_i, xyz_j, xyz_k, xyz_l
   integer(8) i_el, ind_i, ind_j, ind_l, i
   complex(8) p2prod

   double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
   double precision det_b
   double precision p2_p2_p2(6,6,6)
   double precision p2_p2_p2x(6,6,6), p2_p2_p2y(6,6,6)

   bas_ovrlp = C_ZERO

   !       The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
   !       maps the current triangle to the reference triangle.
   do i=1,2
      mat_B(:, i) = nds_xy(:, i+1) - nds_xy(:, 1)
   enddo

   det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)

   if (abs(det_b) .le. 1.0d-22) then
      write(*,*) '?? PE_int_v2: Determinant = 0 :', det_b
      write(*,*) "nds_xy = ", nds_xy
      write(*,*) 'Aborting...'
      stop
   endif

   !       We also need, is the matrix mat_T of the reverse transformation
   !                (from reference to current triangle):
   !       mat_T = inverse matrix of de mat_B
   mat_T(1,1) =  mat_B(2,2) / det_b
   mat_T(2,2) =  mat_B(1,1) / det_b
   mat_T(1,2) = -mat_B(1,2) / det_b
   mat_T(2,1) = -mat_B(2,1) / det_b

   !       Note that if grad_i_0 is the gradient on the reference triangle,
   !       then the gradient on the actual triangle is:
   !       grad_i  = Transpose(mat_T)*grad_i0
   !
   !       mat_T_tr = Transpose(mat_T)
   mat_T_tr(1,1) = mat_T(1,1)
   mat_T_tr(2,2) = mat_T(2,2)
   mat_T_tr(1,2) = mat_T(2,1)
   mat_T_tr(2,1) = mat_T(1,2)

   ! analytic integrals
   call find_overlaps_p2_p2_p2 (p2_p2_p2, det_b)
   call find_overlaps_p2_p2_p2x (p2_p2_p2x, mat_T_tr, det_b)
   call find_overlaps_p2_p2_p2y (p2_p2_p2y, mat_T_tr, det_b)




   do nd_i=1,P2_NODES_PER_EL
      do xyz_i=1,3
         ind_i = xyz_i + 3*(nd_i-1)

         do nd_j=1,P2_NODES_PER_EL
            do xyz_j=1,3
               ind_j = xyz_j + 3*(nd_j-1)

               !  Gradient of transverse components of basis function
               do xyz_k=1,3

                  do nd_l=1,P2_NODES_PER_EL

                     if ( xyz_k .eq. 1) then
                        p2prod = p2_p2_p2x(nd_i,nd_j,nd_l)
                     elseif ( xyz_k .eq. 2) then
                        p2prod = p2_p2_p2y(nd_i,nd_j,nd_l)
                     elseif ( xyz_k .eq. 3) then
                        p2prod = -C_IM_ONE* beta_ac * p2_p2_p2(nd_i,nd_j,nd_l)
                     endif

                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(nd_l-1)
                        bas_ovrlp(ind_i,ind_j,xyz_k,ind_l) = p_tensor(xyz_i,xyz_j,xyz_k,xyz_l,typ_e) * p2prod
                     enddo

                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   bas_ovrlp = bas_ovrlp * eps**2




end subroutine  make_P2_overlaps_i_j_dk_l_analytic


! Calculate the Q_PE integral of two EM modes and an AC mode
! using numerical quadrature over the current triangle
! defined by nds_xy.

! calculates PE style mode overlaps accounting for polarisation
! \int eps_r^2 pe(i,j,k,l) phi2_ir \phi2_js \partial_k \phi_lt  dxdy
! i,j,l run over P2 nodes
! r,s,t run over components x,y,z
!
! k runs over components x,y,z
! \partial k acts as \partial x, \partial y, or \partial z = -i beta
! see storage layout below

subroutine make_P2_overlaps_i_j_dk_l_quadrature(i_el, beta_ac, typ_e, is_curved, &
   eps, p_tensor, n_elt_mats, nds_xy, quadint, bas_ovrlp, errco, emsg)

   use numbatmod
   use class_TriangleIntegrators
   type(QuadIntegrator) quadint

   complex(8) beta_ac, eps
   double precision nds_xy(2,P2_NODES_PER_EL)
   complex(8) p_tensor(3,3,3,3,n_elt_mats)
   integer(8) n_elt_mats, typ_e
   logical is_curved

   ! Contains phi2_i phi2_j \partial_k phi2_l for each basis function
   ! Slots i,j, ,l are indexed phi_1x, phi_1y, phi_1z, phi_2x, phi_2y, ...
   complex(8), dimension(3*P2_NODES_PER_EL, 3*P2_NODES_PER_EL, &
      3, 3*P2_NODES_PER_EL) :: bas_ovrlp

   integer(8), intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   !locals

   double precision t_quadwt
   integer(8) iq, nd_i, nd_j, nd_l, xyz_i, xyz_j, xyz_k, xyz_l
   integer(8) i_el, ind_i, ind_j, ind_l
   double precision phi2_i, phi2_j, phi2_dk_l, phiprod
   complex(8) t_basov
   logical do_P3

   bas_ovrlp = C_ZERO
   ! For each quadrature point evaluate Q_PE of Lagrange polynomials
   ! or derivative of Lagrange polynomials
   do_P3 = .false.

   do iq=1,quadint%n_quad
      call quadint%build_transforms_at(iq, nds_xy, is_curved, do_P3, errco, emsg)
      RETONERROR(errco)

      ! transformed weighting of this quadrature point including triangle area transform
      t_quadwt = quadint%wt_quad(iq) * abs(quadint%det)


      ! Calculate Q_PE of basis functions at quadrature point,
      ! which is a superposition of P2 polynomials for each function (fi_eld).
      do nd_i=1,P2_NODES_PER_EL
         phi2_i = quadint%phi_P2_ref(nd_i)
         do xyz_i=1,3                        ! x/y/z component
            ind_i = xyz_i + 3*(nd_i-1)

            do nd_j=1,P2_NODES_PER_EL
               phi2_j = quadint%phi_P2_ref(nd_j)
               do xyz_j=1,3                      ! x/y/z component
                  ind_j = xyz_j + 3*(nd_j-1)

                  ! Gradient of transverse components of basis function
                  do xyz_k=1,2                    ! x/y component

                     do nd_l=1,P2_NODES_PER_EL
                        phi2_dk_l = quadint%gradt_P2_act(xyz_k,nd_l)
                        phiprod = phi2_i * phi2_j * phi2_dk_l

                        do xyz_l=1,3
                           ind_l = xyz_l + 3*(nd_l-1)

                           t_basov = t_quadwt * p_tensor(xyz_i,xyz_j,xyz_k,xyz_l,typ_e) * phiprod

                           bas_ovrlp(ind_i, ind_j, xyz_k, ind_l) = bas_ovrlp(ind_i, ind_j, xyz_k, ind_l) + t_basov

                        enddo
                     enddo
                  enddo

                  ! Gradient of longitudinal components of basis function,
                  !  which is -i*beta*phi because field is assumed to be of form e^{-i*beta*z} phi.
                  xyz_k=3                        ! z component

                  do nd_l=1,P2_NODES_PER_EL
                     phiprod = phi2_i * phi2_j * quadint%phi_P2_ref(nd_l)

                     do xyz_l=1,3
                        ind_l = xyz_l + 3*(nd_l-1)

                        t_basov = -C_IM_ONE* beta_ac * t_quadwt * p_tensor(xyz_i,xyz_j,xyz_k,xyz_l,typ_e)  * phiprod

                        bas_ovrlp(ind_i,ind_j,xyz_k,ind_l) = bas_ovrlp(ind_i,ind_j,xyz_k,ind_l) + t_basov

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   ! pull the pel_ijkl multiply to here?
   bas_ovrlp = bas_ovrlp * eps**2

end subroutine


