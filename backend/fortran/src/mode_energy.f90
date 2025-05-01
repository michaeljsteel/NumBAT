!
!  Calculate  Int(unit cell) ||E||^2 and Int(over cylinder) epsilon ||E||^2

!  mode_pol(1) contains Int(unit cell) epsilon |E_x|^2 / Int(unit cell) epsilon ||E||^2
!  mode_pol(2) contains Int(unit cell) epsilon |E_y|^2 / Int(unit cell) epsilon ||E||^2
!  mode_pol(3) contains Int(unit cell) epsilon |E_z|^2 / Int(unit cell) epsilon ||E||^2
!  mode_pol(4) contains Int(over cylinder) epsilon ||E||^2 / Int(unit cell) epsilon ||E||^2
!  Note that we use E_z = i * \hat{E}_z * (i*beta) (because of the change of variable)
!  It is assumed that the type number triangle in the cylinders is :
!                           typ_e=n_core(1) or typ_e=n_core(2)
!
subroutine mode_energy (nval, nel,  n_core, &
   mesh_raw, nb_typ_el, eps_eff,  sol, mode_pol)

   use class_MeshRaw
   use numbatmod

   integer(8) nval, nel
   integer(8) nb_typ_el, n_core(2)
   complex(8) sol(3,N_DOF_PER_EL,nval,nel)
   complex(8) eps_eff(nb_typ_el), mode_pol(4,nval)

   type(MeshRawEM) :: mesh_raw

!
!  variables for quadrature interpolation
   integer(8) nquad, nquad_max, iq
   parameter (nquad_max = 25)
   double precision wq(nquad_max)
   double precision xq(nquad_max), yq(nquad_max)
   double precision xx(2), xx_g(2), ww, det
   double precision mat_B(2,2), mat_T(2,2)

   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) vec_phi(3)

   double precision phi2_list(6)

   integer(8) j, iel, ival, typ_e
   integer(8) inode, global, trans
   integer(8) debug, ui

   complex(8) z_tmp, coeff_1


   ui = 6
   debug = 1

   if ( P2_NODES_PER_EL .ne. 6 ) then
      write(ui,*) "overlap_J: problem P2_NODES_PER_EL = ", P2_NODES_PER_EL
      write(ui,*) "overlap_J: P2_NODES_PER_EL should be equal to 6!"
      write(ui,*) "overlap_J: Aborting..."
      stop
   endif

   do ival=1,nval
      do j=1,4
         mode_pol(j,ival) = 0.0d0
      enddo
   enddo

   call quad_triangle(nquad,nquad_max,wq,xq,yq)

!	 loop over all elements

   do iel=1,nel

      typ_e = mesh_raw%el_material(iel)
      do inode=1,P2_NODES_PER_EL
         global = mesh_raw%elnd_to_mshpt(inode,iel)
         nod_el_p(inode) = global
         xel(1,inode) = mesh_raw%v_nd_xy(1,global)
         xel(2,inode) = mesh_raw%v_nd_xy(2,global)
      enddo

      do iq=1,nquad
         xx(1) = xq(iq)
         xx(2) = yq(iq)
         ww = wq(iq)
         call phi2_2d_mat_J(xx, phi2_list)
         call jacobian_p1_2d(xx, xel, P2_NODES_PER_EL, xx_g, det, mat_B, mat_T)

         coeff_1 = ww * ABS(det) * eps_eff(typ_e)
         do ival=1,nval
            do trans=1,3
               vec_phi(trans) = 0.0d0
            enddo
            do inode=1,P2_NODES_PER_EL
               !  transverse field components
               do trans=1,3
                  vec_phi(trans) = vec_phi(trans) + sol(trans,inode,ival,iel) * phi2_list(inode)
               enddo
            enddo

            do trans=1,3
               z_tmp = coeff_1 * abs(vec_phi(trans))**2
               mode_pol(trans,ival) = mode_pol(trans,ival) + z_tmp
            enddo

            if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
               do trans=1,3
                  z_tmp = coeff_1 * abs(vec_phi(trans))**2
                  mode_pol(4,ival) = mode_pol(4,ival) + z_tmp
               enddo
            endif

         enddo
      enddo
   enddo

   !   Total energy and normalization
   do ival=1,nval
      z_tmp = mode_pol(1,ival) + mode_pol(2,ival)  + mode_pol(3,ival)

      !   if (abs(z_tmp) .lt. 1.0d-10) then
      if (abs(z_tmp) .lt. 1.0d-20) then
         write(*,*) "mode_energy: the total energy ",  "is too small : ", z_tmp
         write(*,*) "mode_energy: ival = ", ival
         write(*,*) "mode_energy: zero eigenvector; aborting..."
         stop
      endif

      do j=1,4
         mode_pol(j,ival) = mode_pol(j,ival) / z_tmp
      enddo

   enddo

   return
end
