

!  dof_props = 0  => interiour ddl (ddl = Degree Of Freedom)
!  dof_props != 0 => boundary ddl

!  bdy_cdn = 0 => Dirichlet boundary condition (E-field: electric wall condition)
!  bdy_cdn = 1 => Neumann boundary condition (E-field: magnetic wall condition)
!  bdy_cdn = 2 => Periodic boundary condition



!  This subroutine set the boundary condition parameters


!  The provided Fortran subroutine bound_cond sets boundary condition parameters for a finite element mesh, depending on the specified boundary condition type. Here's a detailed step-by-step explanation of the code:

!  Purpose:
!  The subroutine assigns the degrees of freedom (DOFs) for nodes, edges, and elements in a finite element mesh based on specified boundary conditions (Dirichlet or Neumann).

!  Parameters:
!  Input Parameters:

!    bdy_cdn: Specifies the type of boundary condition.
!      0 for Dirichlet boundary condition (electric wall condition).
!      1 for Neumann boundary condition (magnetic wall condition).
!      2 for Periodic boundary condition.
!
!  n_ddl: Total number of DOFs.

!
!  dof_props: An array where:
!  The first row indicates if a DOF is an interior or boundary DOF.
!  The second row indicates the dimensional type (0 for nodes, 1 for edges, 2 for elements).

!  Output Parameters:

!    n_dof: Total number of equations (DOFs).
!    m_eqs: An array mapping each DOF to its equation number, considering the boundary conditions.

!  Local Variables:
!  i: Loop index.
!  i_boundary: Indicates if the current DOF is on the boundary.
!  i_dim: Dimensional type of the current DOF (node, edge, or element).
!
!  Description:
!    Dirichlet Boundary Condition (bdy_cdn = 0):
!      All points (nodes, edges, elements) have degrees of freedom unless they are boundary points.
!      For elements (i_dim = 2): Each element has 3 interior DOFs.
!      For edges (i_dim = 1): Each edge has 3 DOFs unless it is on the boundary.
!      For nodes (i_dim = 0): Each node has 1 DOF unless it is on the boundary.
!
!  Neumann Boundary Condition (bdy_cdn = 1):
!     All points (nodes, edges, elements) have degrees of freedom.
!     For elements (i_dim = 2) and edges (i_dim = 1): Each has 3 DOFs.
!     For nodes (i_dim = 0): Each node has 1 DOF.

!  matrix m_eqs is built up based on properties of the entities established
!    in MeshEntities
!
! dof_props holds the is_boundary and dimensionality of each entity

!  m_eqs assigns an index to each degree of freedom, if it exists, for each entity

subroutine bound_cond_em (bdy_cdn, n_ddl, dof_props, n_dof, m_eqs, errco, emsg)

   use numbatmod

   integer(8) bdy_cdn, n_ddl, n_dof
   integer(8) m_eqs(3,n_ddl), dof_props(2,n_ddl)
   integer(8),  intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   integer(8) i, is_boundary, i_dim

   if (bdy_cdn .eq. BCS_DIRICHLET) then  !  all points have a degree of freedom

      n_dof = 0
      do i=1,n_ddl
         is_boundary = dof_props(1,i)
         i_dim = dof_props(2,i)

         if (i_dim .eq. 2) then !  each element is associated with 3 interior Degrees Of Freedom (DOF)
            m_eqs(1,i) = n_dof + 1
            m_eqs(2,i) = n_dof + 2
            m_eqs(3,i) = n_dof + 3
            n_dof = n_dof + 3

         elseif (i_dim .eq. 1) then  !  each edge is associated with 3 Degrees Of Freedom (DOF)
            if (is_boundary .eq. 0) then
               m_eqs(1,i) = n_dof + 1
               m_eqs(2,i) = n_dof + 2
               m_eqs(3,i) = n_dof + 3
               n_dof = n_dof + 3
            else                     !  fixed by boundary so no dof
               m_eqs(1,i) = 0
               m_eqs(2,i) = 0
               m_eqs(3,i) = 0
            endif

         elseif (i_dim .eq. 0) then   !  each node is associated with 1 Degree Of Freedom (DOF)
            if (is_boundary .eq. 0) then
               m_eqs(1,i) = n_dof + 1
               m_eqs(2,i) = 0
               m_eqs(3,i) = 0
               n_dof = n_dof + 1
            else
               m_eqs(1,i) = 0
               m_eqs(2,i) = 0
               m_eqs(3,i) = 0
            endif
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
             "bdy_cdn = ", bdy_cdn, "i = ", i
             errco = NBERR_BAD_BOUNDARY_CONDITION
         endif
      enddo

      elseif(bdy_cdn .eq. BCS_NEUMANN) then !all points have a degree of freedom

      n_dof = 0
      do i=1,n_ddl
         i_dim = dof_props(2,i)
         if (i_dim .eq. 2 .or. i_dim .eq. 1) then !  Each element or edge is associated with 3 Degrees Of Freedom (DOF)
            m_eqs(1,i) = n_dof + 1
            m_eqs(2,i) = n_dof + 2
            m_eqs(3,i) = n_dof + 3
            n_dof = n_dof + 3
         elseif (i_dim .eq. 0) then
            m_eqs(1,i) = n_dof + 1
            m_eqs(2,i) = 0
            m_eqs(3,i) = 0
            n_dof = n_dof + 1
         else
            write(emsg,*) "bound_cond: i_dim has invalid value : ", i_dim, &
             "bdy_cdn = ", bdy_cdn, "i = ", i
             errco = NBERR_BAD_BOUNDARY_CONDITION
         endif
      enddo
   endif

end
