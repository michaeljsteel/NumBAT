!  Difference from array_sol_AC.f is that the u_z field is multiplied by i
!  which gives you the correct physical displacement field.

!  v_evecs_arpack(*,i) : contains the imaginary and real parts of the solution for points such that cscmat%m_eqs(i) != 0
!  sol(i) : contains solution for all points
!  The dimension of the geometric domain is : dim_32 = 2
!  The dimension of the vector field is : dim2 = 3


subroutine array_sol_AC (mesh_raw, cscmat, num_modes, iindex, &
 v_eigs_nu, v_evecs_arpack, mode_pol, sol)


   use numbatmod
   use alloc

   use class_MeshRaw
   use class_SparseCSC_AC

   type(MeshRawAC) mesh_raw
   type(SparseCSC_AC) cscmat
   !type(NBError) nberr


   integer(8) num_modes
!  TODO: n_core seems to be never initialised. Is that code ever called?
   integer(8) n_core(2)
   integer(8) iindex(*)
   complex(8) v_evecs_arpack(cscmat%n_dof,num_modes)

!  sol(3, 1..P2_NODES_PER_EL,num_modes, mesh_raw%n_msh_el)          contains the values of the 3 components at P2 interpolation nodes
   complex(8) sol(3,P2_NODES_PER_EL,num_modes,mesh_raw%n_msh_el)
   complex(8) v_eigs_nu(num_modes), v_tmp(num_modes)
   complex(8) mode_pol(4,num_modes)



!  Local variables

   double precision mode_comp(4)
   integer(8) nod_el_p(P2_NODES_PER_EL)
   double precision xel(2,P2_NODES_PER_EL)
   complex(8) sol_el(3,P2_NODES_PER_EL)


   integer(8) j, i1, j1, inod, typ_e, debug
   integer(8) iel, ival, ival2, jp, ind_jp, j_eq
   complex(8) z_tmp1, z_tmp2, z_sol_max

   double precision x_min, x_max, y_min, y_max
   double precision x_mid, y_mid
   double precision dx, dy, x_0, y_0
   double precision lx, ly, rx, ry

   integer(8) i_sol_max, i_sol_max_tmp
   integer(8) i_component, i_component_tmp



   debug = 0


   sol= C_ZERO

   ! do ival=1,num_modes
   !    do iel=1,mesh_raw%n_msh_el
   !       do inod=1,P2_NODES_PER_EL
   !          do j=1,3
   !             sol(j,inod,ival,iel) = 0
   !          enddo
   !       enddo
   !    enddo
   ! enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   x_min = mesh_raw%v_nd_xy(1,1)
   x_max = mesh_raw%v_nd_xy(1,1)
   do j=1,mesh_raw%n_msh_pts
      x_0 = mesh_raw%v_nd_xy(1,j)
      if(x_0 .lt. x_min) x_min = x_0
      if(x_0 .gt. x_max) x_max = x_0
   enddo
   y_min = mesh_raw%v_nd_xy(2,1)
   y_max = mesh_raw%v_nd_xy(2,1)
   do j=1,mesh_raw%n_msh_pts
      y_0 = mesh_raw%v_nd_xy(2,j)
      if(y_0 .lt. y_min) y_min = y_0
      if(y_0 .gt. y_max) y_max = y_0
   enddo

   x_mid = (x_min + x_max) / 2.0d0
   y_mid = (y_min + y_max) / 2.0d0

!  Length in the x direction
   lx = x_max - x_min
!  Length in the y direction
   ly = y_max - y_min

   rx = sqrt(mesh_raw%n_msh_el * lx / (2.0d0 * ly))
   ry = rx * ly/lx
!  rx and ry and such that : rx * ry = 2 * mesh_raw%n_msh_el


   dx = lx / rx
   dy = ly / ry


!c      lx = x_max - x_min  !  Length in the x direction
!c      ly = y_max - y_min  !  Length in the y direction

   if (debug .eq. 1) then
      write(*,*)
      write(*,*) "array_sol: x_min = ", x_min
      write(*,*) "array_sol: x_max = ", x_max
      write(*,*) "array_sol: y_min = ", y_min
      write(*,*) "array_sol: y_max = ", y_max
      write(*,*) "array_sol: x_mid = ", x_mid
      write(*,*) "array_sol: y_mid = ", y_mid
      write(*,*) "array_sol:    dx = ", dx
      write(*,*) "array_sol:    dy = ", dy
      write(*,*) "array_sol:    rx = ", rx
      write(*,*) "array_sol:    ry = ", ry
      write(*,*) "array_sol:   mesh_raw%n_msh_el = ", mesh_raw%n_msh_el
      write(*,*)
   endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   do j=1,num_modes
      j1=iindex(j)
      v_tmp(j) = v_eigs_nu(j1)
   enddo
   do j=1,num_modes
      v_eigs_nu(j) = v_tmp(j)
   enddo

   do ival=1,num_modes
      ival2 = iindex(ival)
      do j=1,4
         mode_pol(j,ival) = 0.0d0
      enddo

      z_sol_max = 0.0d0
      i_sol_max = 0
      i_component = 0

      do iel=1,mesh_raw%n_msh_el
         typ_e = mesh_raw%el_material(iel)
         do j=1,4
            mode_comp(j) = 0.0d0
         enddo
         do inod=1,P2_NODES_PER_EL
            j = mesh_raw%elnd_to_mshpt(inod,iel)
            nod_el_p(inod) = j
            xel(1,inod) = mesh_raw%v_nd_xy(1,j)
            xel(2,inod) = mesh_raw%v_nd_xy(2,j)
         enddo
         do inod=1,P2_NODES_PER_EL
            jp = mesh_raw%elnd_to_mshpt(inod,iel)
            do j_eq=1,3
               ind_jp = cscmat%m_eqs(j_eq,jp)
               if (ind_jp .gt. 0) then
                  z_tmp1 = v_evecs_arpack(ind_jp, ival2)
                  sol_el(j_eq,inod) = z_tmp1
               else
                  sol_el(j_eq,inod) = 0
               endif
            enddo
!  The z-compoenent must be multiplied by ii in order to get the un-normalised z-compoenent
            j_eq=3
            sol_el(j_eq,inod) = C_IM_ONE* sol_el(j_eq,inod)
            do j=1,3
               z_tmp2 = sol_el(j,inod)
               sol(j,inod,ival,iel) = z_tmp2
               if (abs(z_sol_max) .lt. abs(z_tmp2)) then
                  z_sol_max = z_tmp2
!  We want to normalise such the the z-component is purely imaginary complex number
                  if (j == 3) z_sol_max = - C_IM_ONE* z_sol_max
                  i_sol_max = mesh_raw%elnd_to_mshpt(inod,iel)
                  i_component = j
               endif
            enddo
!  Contribution of the element iel to the mode component
            do j=1,3
               z_tmp2 = abs(sol_el(j,inod))**2
               mode_comp(j) = mode_comp(j) + real(z_tmp2)
            enddo
         enddo
!  Avarage values
         do j=1,3
            mode_comp(j) = mode_comp(j)/dble(P2_NODES_PER_EL)
!  mode_comp(j) = abs(det)*mode_comp(j)/dble(P2_NODES_PER_EL)
         enddo
!  Add the contribution of the element iel to the mode component
         do j=1,3
            mode_pol(j,ival) = mode_pol(j,ival) + mode_comp(j)
         enddo

         !  TODO: THIS CODE SEEMS NEVER TO BE CALLED
         !  THESE BAD VALUES TO TRY TO MAKE IT FAIL IF IT DOES
         n_core(1)=-1000
         n_core(2)=-1000

         if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
            write(*,*) "Warning: unitialised n_core in array_sol_AC.f"
            z_tmp2 = mode_comp(1) + mode_comp(2)&
            &+ mode_comp(3)
            mode_pol(4,ival) = mode_pol(4,ival) + z_tmp2
         endif

      enddo
!  Total energy and normalization
      z_tmp2 = mode_pol(1,ival) + mode_pol(2,ival)&
      &+ mode_pol(3,ival)
      if (abs(z_tmp2) .lt. 1.0d-10) then
         write(*,*) "array_sol: the total energy ",&
         &"is too small : ", z_tmp2
         write(*,*) "array_sol: ival ival2 = ", ival, ival2
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif
      do j=1,3
         mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
      enddo
      j=4
      mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
!  Check if the eigenvector is nonzero
      if (abs(z_sol_max) .lt. 1.0d-10) then
         z_sol_max = z_tmp2
         write(*,*) "array_sol: z_sol_max is too small"
         write(*,*) "array_sol: z_sol_max = ", z_sol_max
         write(*,*) "ival, ival2, num_modes = ", ival, ival2, num_modes
         write(*,*) "array_sol: zero eigenvector; aborting..."
         stop
      endif
      if (debug .eq. 1) then
         write(*,*) "array_sol (A): "
         write(*,*) "                           ival = ", ival
         write(*,*) "z_sol_max     = ", z_sol_max, abs(z_sol_max)
      endif
!  Normalization so that the maximum field component is 1
      do iel=1,mesh_raw%n_msh_el
         do inod=1,P2_NODES_PER_EL
            i1 = mesh_raw%elnd_to_mshpt(inod,iel)
            do j=1,3
               z_tmp1 = sol(j,inod,ival,iel)/z_sol_max
               sol(j,inod,ival,iel) = z_tmp1
            enddo
            i1 = mesh_raw%elnd_to_mshpt(inod,iel)
            if (i1 .eq. i_sol_max_tmp .and. debug .eq. 1) then
               write(*,*) "array_sol (B):"
               write(*,*) "ival, i1, iel = ", ival, i1, iel
               write(*,*) "array_sol: Field normalisaion point:"
               write(*,*) "x = ", dble(mesh_raw%v_nd_xy(1,i1))
               write(*,*) "y = ", dble(mesh_raw%v_nd_xy(2,i1))
               write(*,*) "i_sol_max = ", i_sol_max
               write(*,*) "i_component_tmp = ", i_component_tmp
               write(*,*) ival, i1, iel,&
               &(dble(sol(j,inod,ival,iel)),j=1,3)
               write(*,*) ival, i1, iel,&
               &(imag(sol(j,inod,ival,iel)),j=1,3)
               write(*,*) ival, i1, iel,&
               &(abs(sol(j,inod,ival,iel)),j=1,3)
            endif
         enddo
      enddo
      do j=1,cscmat%n_dof
         z_tmp1 = v_evecs_arpack(j,ival2)/z_sol_max
         v_evecs_arpack(j,ival2) = z_tmp1
      enddo

      if (debug .eq. 1) then
         write(*,*)
      endif
   enddo

   if (debug .eq. 1) then
      write(*,*) "Exiting subroutine array_sol"
   endif



   return
end
