#include "numbat_decl.h"

module class_MeshProps

   use numbatmod
   use alloc

   implicit none


   private

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

      !character(FNAME_LENGTH) :: mesh_file

      integer(8) :: n_msh_pts
      integer(8) :: n_msh_el
      integer(8) :: n_elt_mats
      integer(8) :: n_ddl  ! = 9 * n_msh_el

      integer(8) :: nodes_per_el


      integer(8), dimension(:), allocatable :: type_nod          ! rename to node_phys_i
      double precision, dimension(:,:), allocatable :: xy_nodes  ! rename to

      integer(8), dimension(:), allocatable :: type_el           ! rename to el_mat_i
      integer(8), dimension(:,:), allocatable :: table_nod       ! .mail elt data, but transposed: tn[i,j] = ith node of elt j
      !  tn[1..3,j] = vertex nodes, [4..6, j] edge nodes

   contains

      procedure :: init => MeshProps_init
      !final :: MeshProps_destructor

      procedure :: fill_python_arrays => MeshProps_fill_python_arrays
      procedure :: is_boundary_node => MeshProps_is_boundary_node
      procedure :: is_boundary_node_2 => MeshProps_is_boundary_node_2

      procedure :: type_node_by_ref => MeshProps_type_node_by_ref



   end type MeshProps

! ---------------------------------------


   type, public  :: N_E_F_Props


      ! table_nod[i, j] 14 x n_msh_el
      !  tn[1, j] = j  ! index of element
      !    -


      integer(8), dimension(:,:), allocatable :: type_nod
      integer(8), dimension(:,:), allocatable :: table_nod
      double precision, dimension(:,:), allocatable :: xy_nodes

   contains

      procedure :: init => N_E_F_Props_init
      ! final :: N_E_F_Props_destructor

   end type N_E_F_Props

   ! ---------------------------------------
contains

   subroutine MeshProps_init(this, n_msh_pts, n_msh_el, n_elt_mats, &
      errco, emsg)

      class(MeshProps) :: this
      integer(8) :: n_msh_el, n_msh_pts, n_elt_mats
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      this%n_msh_pts = n_msh_pts
      this%n_msh_el = n_msh_el
      this%n_ddl = 9 * n_msh_el
      this%n_elt_mats = n_elt_mats

      this%nodes_per_el = 6


      call integer_alloc_1d(this%type_el, n_msh_el, 'type_el', errco, emsg); RETONERROR(errco)

      call double_alloc_2d(this%xy_nodes, 2_8, n_msh_pts, 'xy_nodes', errco, emsg); RETONERROR(errco)

      call integer_alloc_1d(this%type_nod, n_msh_pts, 'type_nod', errco, emsg); RETONERROR(errco)

      call integer_alloc_2d(this%table_nod, this%nodes_per_el, n_msh_el, 'table_nod', errco, emsg);
      RETONERROR(errco)

   end subroutine

   ! subroutine MeshProps_destructor(this)
   !    type(MeshProps) :: this

   ! end subroutine

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

   ! boundary nodes have non zero GMsh physindex codes
   pure logical function  MeshProps_is_boundary_node(this, nd) result(res)
      class(MeshProps), intent(in) :: this
      integer(8), intent(in)  :: nd

      res = this%type_nod(nd) .ne. 0
   end function

   pure logical function MeshProps_is_boundary_node_2(this, i_nd, i_el) result(res)
      class(MeshProps), intent(in) :: this
      integer(8), intent(in)  :: i_nd, i_el
      res = this%type_nod(this%table_nod(i_nd, i_el)) .ne. 0
   end function

   ! get node type by indirection through node table
   integer(8) function  MeshProps_type_node_by_ref(this, i_nd, i_el) result(res)
      class(MeshProps) :: this
      integer(8) :: i_nd, i_el
      res = this%type_nod(this%table_nod(i_nd, i_el))
   end function






   subroutine N_E_F_Props_init(this, n_msh_el, n_ddl, &
      errco, emsg)

      class(N_E_F_Props) :: this
      integer(8) :: n_msh_el, n_ddl
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      call double_alloc_2d(this%xy_nodes, 2_8, n_ddl, 'xy_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(this%type_nod, 2_8, n_ddl, 'type_N_E_F', errco, emsg); RETONERROR(errco)
      call integer_alloc_2d(this%table_nod, 14_8, n_msh_el, 'table_N_E_F', errco, emsg); RETONERROR(errco)

   end subroutine

   ! subroutine N_E_F_Props_destructor(this)
   !    type(N_E_F_Props) :: this

   ! end subroutine




end module class_MeshProps
