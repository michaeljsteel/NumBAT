#include "numbat_decl.h"


 ! Calculate the H-field soln_H1 from the E-field soln_k1 of a mode
 !  The z-component of the E-field is not normalised

subroutine H_mode_field_Ez (k_0, n_modes, n_msh_el, n_msh_pts, nnodes_P2, elnd_to_mshpt, &
 x, v_beta, soln_k1, soln_H1)

   use numbatmod
   integer(8) n_modes, n_msh_el, n_msh_pts, nnodes_P2
   integer(8) elnd_to_mshpt(nnodes_P2,n_msh_el)
   double precision x(2,n_msh_pts)
   complex(8) soln_k1(3,nnodes_P2+7,n_modes,n_msh_el)
   complex(8) soln_H1(3,nnodes_P2,n_modes,n_msh_el)
   complex(8) beta1
   complex(8) v_beta(n_modes)
   double precision k_0

   !     Local variables

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) E_field_el(3,P2_NODES_PER_EL)
   complex(8) H_field_el(3,P2_NODES_PER_EL)

   !   P3 Ez-field
   complex(8) Ez_field_el_P3(P3_NODES_PER_EL)
   integer(8) i, j, j1
   integer(8) iel, ival, inod
   integer(8) ui
   !complex(8) z_tmp1

   double precision mat_B(2,2), mat_T(2,2), det_b
   integer(8), parameter :: ZCOMP = 3

   !f2py intent(in) k_0, n_modes, n_msh_el, n_msh_pts
   !f2py intent(in) nnodes_P2, elnd_to_mshpt
   !f2py intent(in) x, v_beta, soln_k1
   !
   !f2py depend(elnd_to_mshpt) nnodes_P2, n_msh_el
   !f2py depend(x) n_msh_pts
   !f2py depend(v_beta) n_modes
   !f2py depend(soln_k1) nnodes_P2, n_modes, n_msh_el
   !
   !f2py intent(out) soln_H1

   ui = stdout


   if ( nnodes_P2 .ne. 6 ) then
      write(ui,*) "EM_mode_en_int_v2: problem nnodes = ", nnodes_P2
      write(ui,*) "EM_mode_en_int_v2: nnodes should be equal to 6 !"
      write(ui,*) "EM_mode_en_int_v2: Aborting..."
      stop
   endif

   do ival=1,n_modes
      beta1 = v_beta(ival)
      do iel=1,n_msh_el
         do j=1,nnodes_P2
            j1 = elnd_to_mshpt(j,iel)
            nod_el_p(j) = j1
            xel(1,j) = x(1,j1)
            xel(2,j) = x(2,j1)
         enddo

         !         The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
         !         maps the current triangle to the reference triangle.
         do i=1,2
            do j=1,2
               mat_B(j,i) = xel(j,i+1) - xel(j,1)
            enddo
         enddo

         det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
         if (abs(det_b) .le. 1.0d-22) then
            write(*,*) '?? H_mode_field_Ez: Deter. = 0 :', det_b
            write(*,*) "xel = ", xel
            write(*,*) 'Aborting...'
            stop
         endif

         !         We also need  the matrix mat_T of the reverse transmation
         !                  (from reference to current triangle):
         !         mat_T = inverse matrix of de mat_B
         mat_T(1,1) =  mat_B(2,2) / det_b
         mat_T(2,2) =  mat_B(1,1) / det_b
         mat_T(1,2) = -mat_B(1,2) / det_b
         mat_T(2,1) = -mat_B(2,1) / det_b
         !         Note that if grad_i_0 is the gradient on the reference triangle,
         !         then the gradient on the actual triangle is:
         !         grad_i  = Transpose(mat_T)*grad_i0

         do inod=1,nnodes_P2
            !           The components (E_x,E_y) of the mode ival
            do j=1,2

               E_field_el(j,inod) = soln_k1(j,inod,ival,iel)
            enddo

            !           The component E_z of the mode ival. The FEM code uses the scaling:
            !           E_z = C_IM_ONE* beta1 * \hat{E}_z
            !j=3

            E_field_el(j,inod) = soln_k1(ZCOMP,inod,ival,iel)
         enddo

         !         E_z-field:
         do inod=1,3
            !           The longitudinal component at the vertices (P3 elements)

            Ez_field_el_P3(inod) = soln_k1(ZCOMP,inod,ival,iel)
         enddo

         do inod=4,P3_NODES_PER_EL
            !           The longitudinal component at the edge nodes and interior node (P3 elements)
            !j=3

            Ez_field_el_P3(inod) = soln_k1(ZCOMP,inod+nnodes_P2-3,ival,iel)
         enddo

         call get_H_field_p3(k_0, beta1, mat_T,&
          E_field_el, Ez_field_el_P3, H_field_el)

         do inod=1,P2_NODES_PER_EL
            do j=1,3
               soln_H1(j,inod,ival,iel) = H_field_el(j,inod)
            enddo
         enddo
      enddo
   enddo

end subroutine H_mode_field_Ez
