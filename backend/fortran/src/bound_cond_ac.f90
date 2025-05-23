

!  mesh%v_mshpt_physindex = 0  => interior node
!  mesh%v_mshpt_physindex != 0 => boundary node

!  bdy_cdn = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!  bdy_cdn = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!  bdy_cdn = 2 => Periodic boundary condition



!  This subroutine set the boundary condition parameters
!subroutine bound_cond_AC (bdy_cdn, npt, n_dof, type_nod, in_dof)

subroutine bound_cond_AC (bdy_cdn, mesh, n_dof, in_dof)

   use numbatmod
   use class_Mesh
   type(MeshAC) mesh

   integer(8) bdy_cdn, n_dof
   integer(8) in_dof(3,mesh%n_msh_pts)
   integer(8) i
   logical is_interior

   if(bdy_cdn .eq. BCS_DIRICHLET) then   !  all interior points have a degree of freedom

       n_dof = 0
      do i=1,mesh%n_msh_pts
         is_interior = mesh%v_mshpt_physindex(i) == 0

         if (is_interior ) then !  each element is associated to 3 interior DOF
            in_dof(1,i) = n_dof + 1
            in_dof(2,i) = n_dof + 2
            in_dof(3,i) = n_dof + 3
            n_dof = n_dof + 3
         else
            in_dof(1,i) = 0
            in_dof(2,i) = 0
            in_dof(3,i) = 0
         endif
      enddo
   elseif(bdy_cdn .eq. BCS_NEUMANN) then !  all points have a degree of freedom
      n_dof = 0
      do i=1,mesh%n_msh_pts
         in_dof(1,i) = n_dof + 1
         in_dof(2,i) = n_dof + 2
         in_dof(3,i) = n_dof + 3
         n_dof = n_dof + 3
      enddo
   endif

   return
end
