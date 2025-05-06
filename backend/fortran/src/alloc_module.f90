#include "numbat_decl.h"

module alloc
   use numbatmod

contains


   subroutine logical_alloc_1d(vec, m, nm, nberr)
      logical, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer4_alloc_1d(vec, m, nm, nberr)
      integer(4), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer_alloc_1d(vec, m, nm, nberr)
      integer(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine


   subroutine double_alloc_1d(vec, m, nm, nberr)
      double precision, dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)

      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)

   end subroutine

   subroutine complex_alloc_1d(vec, m, nm, nberr)
      complex(8), dimension(:), allocatable, intent(inout) :: vec
      integer(8) m
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m), STAT=alloc_stat)
      call check_alloc(alloc_stat, m, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer_alloc_2d(vec, m, n, nm, nberr)
      integer(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine integer4_alloc_2d(vec, m, n, nm, nberr)
      integer(4), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m, n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine double_alloc_2d(vec, m, n, nm, nberr)
      double precision, dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine complex_alloc_2d(vec, m, n, nm, nberr)
      complex(8), dimension(:,:), allocatable, intent(inout) :: vec
      integer(8) m,n
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n, nm, NBERROR_130, nberr)
   end subroutine

   subroutine complex_alloc_4d(vec, m, n, q, r, nm, nberr)
      complex(8), dimension(:,:,:,:), allocatable, intent(inout) :: vec
      integer(8) m, n, q, r
      character(*) nm

      integer(8) alloc_stat
      type(NBError) nberr

      call nberr%reset()

      allocate(vec(m,n,q,r), STAT=alloc_stat)
      call check_alloc(alloc_stat, m*n*q*r, nm, NBERROR_130, nberr)
   end subroutine






end module

