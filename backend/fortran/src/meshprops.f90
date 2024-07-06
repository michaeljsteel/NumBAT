#include "numbat_decl.h"

module class_MeshProps

   use numbatmod
   use alloc

   implicit none
   private

   integer(8), parameter :: nodes_per_el = 6

   ! ---------------------------------------

   type, public  :: MeshCounts


      character(FNAME_LENGTH) :: mesh_file
      integer(8) :: n_msh_pts
      integer(8) :: n_msh_el
      integer(8) :: n_elt_mats
      integer(8) :: n_ddl  ! = 9 * n_msh_el


   end type MeshCounts

   ! ---------------------------------------

   type, public  :: MeshProps

   integer(8), dimension(:), allocatable :: type_el
   integer(8), dimension(:), allocatable :: type_nod
   integer(8), dimension(:,:), allocatable :: table_nod
   double precision, dimension(:,:), allocatable :: xy_nodes

contains

   procedure :: init => MeshProps_init
   final :: MeshProps_destructor

   procedure :: fill_python_arrays => MeshProps_fill_python_arrays

end type MeshProps

! ---------------------------------------


   type, public  :: N_E_F_Props


      integer(8), dimension(:,:), allocatable :: type_nod
      integer(8), dimension(:,:), allocatable :: table_nod
      double precision, dimension(:,:), allocatable :: xy_nodes

   contains

      procedure :: init => N_E_F_Props_init
      final :: N_E_F_Props_destructor

   end type N_E_F_Props

   ! ---------------------------------------
contains

subroutine MeshProps_init(this, n_msh_pts, n_msh_el, &
    errco, emsg)

      class(MeshProps) :: this
      integer(8) :: n_msh_el, n_msh_pts
      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      call integer_alloc_1d(this%type_el, n_msh_el, 'type_el', errco, emsg); RETONERROR(errco)

      call double_alloc_2d(this%xy_nodes, 2_8, n_msh_pts, 'xy_nodes', errco, emsg); RETONERROR(errco)

      call integer_alloc_1d(this%type_nod, n_msh_pts, 'type_nod', errco, emsg); RETONERROR(errco)

      call integer_alloc_2d(this%table_nod, nodes_per_el, n_msh_el, 'table_nod', errco, emsg); RETONERROR(errco)

   end subroutine

   subroutine MeshProps_destructor(this)
      type(MeshProps) :: this

   end subroutine

   subroutine MeshProps_fill_python_arrays(this, &
    type_el, type_nod, table_nod, xy_nodes)

    class(MeshProps) :: this

    integer(8), intent(out) :: type_el(:)
    integer(8), intent(out) :: type_nod(:)
    integer(8), intent(out) :: table_nod(:, :)
    double precision, intent(out) :: xy_nodes(:,:)

    type_el = this%type_el
    type_nod = this%type_nod
    table_nod = this%table_nod
    xy_nodes = this%xy_nodes

   end subroutine



subroutine N_E_F_Props_init(this, n_msh_el, n_ddl, &
    errco, emsg)

      class(N_E_F_Props) :: this
      integer(8) :: n_msh_el, n_ddl
      integer, intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      call double_alloc_2d(this%xy_nodes, 2_8, n_ddl, 'xy_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(this%type_nod, 2_8, n_ddl, 'type_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(this%table_nod, 14_8, n_msh_el, 'table_N_E_F', errco, emsg); RETONERROR(errco)

   end subroutine

   subroutine N_E_F_Props_destructor(this)
      type(N_E_F_Props) :: this

   end subroutine




end module class_MeshProps
