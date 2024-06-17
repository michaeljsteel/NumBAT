#include "numbat_decl.h"

subroutine prepare_workspaces(is_em, n_msh_pts, n_msh_el, n_modes, &
   int_max, cmplx_max, real_max, &
   a_iwork, b_zwork, c_dwork, d_dwork, iindex, overlap_L, &
   errco, emsg)

   use numbatmod

   integer :: is_em
   integer*8 :: int_max, cmplx_max, real_max
   integer*8 :: n_msh_el, n_msh_pts
   integer*8 :: n_modes

   integer*8, dimension(:), allocatable :: a_iwork
   complex*16, dimension(:), allocatable :: b_zwork
   double precision, dimension(:), allocatable :: c_dwork
   double precision, dimension(:,:), allocatable :: d_dwork

   integer*8, dimension(:), allocatable :: iindex
   complex*16, dimension(:,:), allocatable :: overlap_L
   integer :: errco
   character(len=EMSG_LENGTH) :: emsg

   integer*8 n_ddl
   integer :: alloc_stat=0

   call array_size(n_msh_pts, n_msh_el, n_modes, &
      int_max, cmplx_max, real_max, n_ddl, errco, emsg)
   RETONERROR(errco)

   allocate(a_iwork(int_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, int_max, "a", -1, errco, emsg)
   RETONERROR(errco)

   allocate(b_zwork(cmplx_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, cmplx_max, "b", -1, errco, emsg)
   RETONERROR(errco)

   allocate(c_dwork(real_max), STAT=alloc_stat)
   call check_alloc(alloc_stat, real_max, "c", -1, errco, emsg)
   RETONERROR(errco)

   allocate(iindex(n_modes), STAT=alloc_stat)
   call check_alloc(alloc_stat, n_modes, "iindex", -1, errco, emsg)
   RETONERROR(errco)

   if (is_em > 0) then
      allocate(d_dwork(2,n_ddl), STAT=alloc_stat)
      call check_alloc(alloc_stat, 2*n_ddl, "d_dwork", -1, errco, emsg)
      RETONERROR(errco)

      allocate(overlap_L(n_modes,n_modes), STAT=alloc_stat)
      call check_alloc(alloc_stat, n_modes*n_modes, "overlap_L", -1, errco, emsg)
      RETONERROR(errco)


   endif

end subroutine prepare_workspaces

!   ----------------------------------------------------------------------------------------

subroutine set_boundary_conditions(bdy_cdn, n_msh_pts, n_msh_el, mesh_xy, nodes_per_el, &
   type_nod, table_nod, n_ddl, neq, ip_type_N_E_F, ip_eq, &
   a_iwork, d_dwork, int_max, debug)

   use numbatmod

   integer*8 :: bdy_cdn, neq, n_msh_pts, n_msh_el, n_ddl, nodes_per_el
   integer*8 :: ip_type_N_E_F, ip_eq, jp_x_n_e_f, int_max
   integer*8 :: ip_period_N, ip_nperiod_N, ip_period_N_E_F, ip_nperiod_N_E_F
   integer*8 :: debug
   double precision mesh_xy(2,n_msh_pts)
   integer*8 type_nod(n_msh_pts)
   integer*8 table_nod(nodes_per_el, n_msh_el)

   ! is this the right way to pass these?
   integer*8, dimension(int_max) :: a_iwork
   double precision, dimension(2,n_ddl) :: d_dwork

   double precision, dimension(2,2) :: lat_vecs

   if ( bdy_cdn .eq. BCS_DIRICHLET .or.  bdy_cdn .eq. BCS_NEUMANN) then

      call bound_cond ( bdy_cdn, n_ddl, neq, a_iwork(ip_type_N_E_F), a_iwork(ip_eq))

   elseif( bdy_cdn .eq. BCS_PERIODIC) then  !  Periodic  conditions (never in NumBAT)
      if (debug .eq. 1) then
         write(*,*) "###### periodic_node"
      endif

      ! reproduced in py_calc_modes.f
      jp_x_N_E_F = 1
      ip_period_N = ip_type_N_E_F + 2*n_ddl
      ip_nperiod_N = ip_period_N + n_msh_pts
      ip_period_N_E_F = ip_nperiod_N + n_msh_pts
      ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
      ip_eq = ip_nperiod_N_E_F + n_ddl

      call lattice_vec (n_msh_pts, mesh_xy, lat_vecs, debug)

      call periodic_node(n_msh_el, n_msh_pts, nodes_per_el, type_nod, mesh_xy, a_iwork(ip_period_N), &
         a_iwork(ip_nperiod_N), table_nod, lat_vecs)

      if (debug .eq. 1) then
         write(*,*) "set_boundary_conditions: ###### periodic_N_E_F"
      endif

      call periodic_N_E_F (n_ddl, a_iwork(ip_type_N_E_F), d_dwork, a_iwork(ip_period_N_E_F), &
         a_iwork(ip_nperiod_N_E_F), lat_vecs)

      call periodic_cond ( bdy_cdn, n_ddl, neq, a_iwork(ip_type_N_E_F), &
         a_iwork(ip_period_N_E_F), a_iwork(ip_eq), debug)

      if (debug .eq. 1) then
         write(*,*) "py_calc_modes.f: neq, n_ddl = ", neq, n_ddl
      endif
      write(*,*) "py_calc_modes.f: neq, n_ddl = ", neq, n_ddl

   endif

end subroutine set_boundary_conditions
