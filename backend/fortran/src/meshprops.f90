#include "numbat_decl.h"

module class_MeshProps

   use numbatmod
   use alloc

   implicit none
   private

   integer(8), parameter :: nodes_per_el = 6

   ! ---------------------------------------

   type, public  :: MeshProps


      character(FNAME_LENGTH) :: mesh_file
      integer(8) :: n_msh_pts
      integer(8) :: n_msh_el
      integer(8) :: n_elt_mats

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
