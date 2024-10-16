

!  type_nod = 0  => interior node
!  type_nod != 0 => boundary node

!  i_cond = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!  i_cond = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!  i_cond = 2 => Periodic boundary condition



!  This subroutine set the boundary condition parameters

subroutine bound_cond_AC (i_cond, npt, n_dof, type_nod, in_dof)

   use numbatmod
   integer(8) i_cond, npt, n_dof
   integer(8) in_dof(3,npt), type_nod(npt)

   integer(8) i
   logical is_interior

   if(i_cond .eq. BCS_DIRICHLET) then   !  all interior points have a degree of freedom

       n_dof = 0
      do i=1,npt
         is_interior = type_nod(i) == 0

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
   elseif(i_cond .eq. BCS_NEUMANN) then !  all points have a degree of freedom
      n_dof = 0
      do i=1,npt
         in_dof(1,i) = n_dof + 1
         in_dof(2,i) = n_dof + 2
         in_dof(3,i) = n_dof + 3
         n_dof = n_dof + 3
      enddo
   endif

   return
end
