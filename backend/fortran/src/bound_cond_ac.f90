
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

subroutine bound_cond_AC (i_cond, npt, neq, type_nod, ineq)

   use numbatmod
   integer(8) i_cond, npt, neq
   integer(8) ineq(3,npt), type_nod(npt)

   integer(8) i
   logical is_interior

   if(i_cond .eq. BCS_DIRICHLET) then   ! all interior points have a degree of freedom

       neq = 0
      do i=1,npt
         is_interior = type_nod(i) == 0

         if (is_interior ) then ! each element is associated to 3 interior DOF
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
   elseif(i_cond .eq. BCS_NEUMANN) then ! all points have a degree of freedom
      neq = 0
      do i=1,npt
         ineq(1,i) = neq + 1
         ineq(2,i) = neq + 2
         ineq(3,i) = neq + 3
         neq = neq + 3
      enddo
   endif
!
   return
end
