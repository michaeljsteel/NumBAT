#include "numbat_decl.h"

module alloc
   use numbatmod

contains

   subroutine logical_alloc_1d(vec, m, nm, errco, emsg)
      logical, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine integer4_alloc_1d(vec, m, nm, errco, emsg)
      integer(4), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine integer_alloc_1d(vec, m, nm, errco, emsg)
      integer(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine


   subroutine double_alloc_1d(vec, m, nm, errco, emsg)
      double precision, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)

      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)

   end subroutine

   subroutine complex_alloc_1d(vec, m, nm, errco, emsg)
      complex(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine integer_alloc_2d(vec, m, n, nm, errco, emsg)
      integer(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine integer4_alloc_2d(vec, m, n, nm, errco, emsg)
      integer(4), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine double_alloc_2d(vec, m, n, nm, errco, emsg)
      double precision, dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine complex_alloc_2d(vec, m, n, nm, errco, emsg)
      complex(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine

   subroutine complex_alloc_4d(vec, m, n, q, r, nm, errco, emsg)
      complex(8), dimension(:,:,:,:), allocatable, intent(inout) :: vec
      integer(8) m, n, q, r
      character(*) nm
      integer(8),  intent(out) :: errco
      character(len=EMSG_LENGTH), intent(out) :: emsg

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n,q,r), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n*q*r, nm, NBERROR_130, nberr)
      call nberr%to_py(errco, emsg)
   end subroutine











   subroutine logical_nalloc_1d(vec, m, nm, nberr)
      logical, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer4_nalloc_1d(vec, m, nm, nberr)
      integer(4), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer_nalloc_1d(vec, m, nm, nberr)
      integer(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine


   subroutine double_nalloc_1d(vec, m, nm, nberr)
      double precision, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m), STAT=alloc_stat)

      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)

   end subroutine

   subroutine complex_nalloc_1d(vec, m, nm, nberr)
      complex(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer_nalloc_2d(vec, m, n, nm, nberr)
      integer(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer4_nalloc_2d(vec, m, n, nm, nberr)
      integer(4), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine double_nalloc_2d(vec, m, n, nm, nberr)
      double precision, dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine complex_nalloc_2d(vec, m, n, nm, nberr)
      complex(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine complex_nalloc_4d(vec, m, n, q, r, nm, nberr)
      complex(8), dimension(:,:,:,:), allocatable, intent(inout) :: vec
      integer(8) m, n, q, r
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      allocate(vec(m,n,q,r), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n*q*r, nm, NBERROR_130, nberr)
   end subroutine






end module

