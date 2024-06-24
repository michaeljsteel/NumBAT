
!
!   type_nod = 0  => interior node
!   type_nod != 0 => boundary node
!
!   i_cond = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!   i_cond = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!   i_cond = 2 => Periodic boundary condition
!
!

!     This subroutine set the boundary condition parameters

subroutine bound_cond_AC (i_cond, npt, neq, type_nod, ineq, debug)

   use numbatmod

   integer(8) i_cond, npt, neq, debug
   integer(8) ineq(3,npt), type_nod(npt)

   integer(8) i, i_boundary

   if(i_cond .eq. BCS_DIRICHLET) then

      !       Dirichlet boundary condition (fixed interface): all interior points have a degree of freedom
      neq = 0
      do i=1,npt
         i_boundary = type_nod(i)
!           ! each element is associated to 3 interior Degrees Of Freedom (DOF)
         if (i_boundary .eq. 0) then
            ineq(1,i) = neq + 1
            ineq(2,i) = neq + 2
            ineq(3,i) = neq + 3
            neq = neq + 3
         else
            ineq(1,i) = 0
            ineq(2,i) = 0
            ineq(3,i) = 0
         endif
      enddo

   elseif(i_cond .eq. BCS_NEUMANN) then
!       Neumann boundary condition (free interface): all points have a degree of freedom

      neq = 0
      do i=1,npt
         ineq(1,i) = neq + 1
         ineq(2,i) = neq + 2
         ineq(3,i) = neq + 3
         neq = neq + 3
      enddo
   endif

   return
end
