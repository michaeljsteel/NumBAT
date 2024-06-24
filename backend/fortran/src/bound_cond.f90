
!
!   type_N_E_F = 0  => interiour ddl (ddl = Degree Of Freedom)
!   type_N_E_F != 0 => boundary ddl
!
!   bdy_cdn = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!   bdy_cdn = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!   bdy_cdn = 2 => Periodic boundary condition
!
!

!     This subroutine set the boundary condition parameters


! The provided Fortran subroutine bound_cond sets boundary condition parameters for a finite element mesh, depending on the specified boundary condition type. Here's a detailed step-by-step explanation of the code:

! Purpose:
! The subroutine assigns the degrees of freedom (DOFs) for nodes, edges, and elements in a finite element mesh based on specified boundary conditions (Dirichlet or Neumann).

! Parameters:
! Input Parameters:

! bdy_cdn: Specifies the type of boundary condition.
! 0 for Dirichlet boundary condition (electric wall condition).
! 1 for Neumann boundary condition (magnetic wall condition).
! 2 for Periodic boundary condition.
! n_ddl: Total number of DOFs.
! type_N_E_F: An array where:
! The first row indicates if a DOF is an interior or boundary DOF.
! The second row indicates the dimensional type (0 for nodes, 1 for edges, 2 for elements).
! Output Parameters:

! neq: Total number of equations (DOFs).
! ineq: An array mapping each DOF to its equation number, considering the boundary conditions.
! Local Variables:
! i: Loop index.
! i_boundary: Indicates if the current DOF is on the boundary.
! i_dim: Dimensional type of the current DOF (node, edge, or element).
! Description:
! Dirichlet Boundary Condition (bdy_cdn = 0):

! All points (nodes, edges, elements) have degrees of freedom unless they are boundary points.
! For elements (i_dim = 2): Each element has 3 interior DOFs.
! For edges (i_dim = 1): Each edge has 3 DOFs unless it is on the boundary.
! For nodes (i_dim = 0): Each node has 1 DOF unless it is on the boundary.
! Neumann Boundary Condition (bdy_cdn = 1):

! All points (nodes, edges, elements) have degrees of freedom.
! For elements (i_dim = 2) and edges (i_dim = 1): Each has 3 DOFs.
! For nodes (i_dim = 0): Each node has 1 DOF.
! Error Handling:

! If bdy_cdn has an invalid value, or if i_dim is not 0, 1, or 2, the subroutine prints an error message and stops execution.



subroutine bound_cond (bdy_cdn, n_ddl, neq, type_N_E_F, ineq)

   use numbatmod

   integer(8) bdy_cdn, n_ddl, neq
   integer(8) ineq(3,n_ddl), type_N_E_F(2,n_ddl)

   integer(8) i, i_boundary, i_dim

   if(bdy_cdn .eq. BCS_DIRICHLET) then  ! all points have a degree of freedom

      neq = 0
      do i=1,n_ddl
         i_boundary = type_N_E_F(1,i)
         i_dim = type_N_E_F(2,i)

         if (i_dim .eq. 2) then ! each element is associated to 3 interior Degrees Of Freedom (DOF)
            ineq(1,i) = neq + 1
            ineq(2,i) = neq + 2
            ineq(3,i) = neq + 3
            neq = neq + 3

         elseif (i_dim .eq. 1) then  ! each edge is associated to 3 Degrees Of Freedom (DOF)
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

         elseif (i_dim .eq. 0) then   ! each nodee is associated to 1 Degree Of Freedom (DOF)
            if (i_boundary .eq. 0) then
               ineq(1,i) = neq + 1
               ineq(2,i) = 0
               ineq(3,i) = 0
               neq = neq + 1
            else
               ineq(1,i) = 0
               ineq(2,i) = 0
               ineq(3,i) = 0
            endif
         else
            write(*,*) "bound_cond: i_dim has invalid value : ", i_dim
            write(*,*) "bound_cond: bdy_cdn = ", bdy_cdn
            write(*,*) "bound_cond: i = ", i
            write(*,*) "bound_cond: Aborting..."
            stop
         endif
      enddo

      elseif(bdy_cdn .eq. BCS_NEUMANN) then !all points have a degree of freedom

      neq = 0
      do i=1,n_ddl
         i_dim = type_N_E_F(2,i)
         if (i_dim .eq. 2 .or. i_dim .eq. 1) then ! Each element or edge is associated to 3 Degrees Of Freedom (DOF)
            ineq(1,i) = neq + 1
            ineq(2,i) = neq + 2
            ineq(3,i) = neq + 3
            neq = neq + 3
         elseif (i_dim .eq. 0) then
            ineq(1,i) = neq + 1
            ineq(2,i) = 0
            ineq(3,i) = 0
            neq = neq + 1
         else
            write(*,*) "bound_cond: i_dim has invalid value : ", i_dim
            write(*,*) "bound_cond: bdy_cdn = ", bdy_cdn
            write(*,*) "bound_cond: i = ", i
            write(*,*) "bound_cond: Aborting..."
            stop
         endif
      enddo
   else
      write(*,*) "bound_cond: bdy_cdn has invalid value : ", bdy_cdn
      write(*,*) "bound_cond: Aborting..."
      stop
   endif

   return
end
