#include "numbat_decl.h"

module alloc
    use numbatmod

contains

    subroutine logical_alloc_1d(vec, m, nm, errco, emsg)
        logical, dimension(:), allocatable, intent(inout) :: vec
        integer(8) m
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m), STAT=alloc_stat)
        call check_alloc(alloc_stat, m, nm, NBERROR_130, errco, emsg)
        end subroutine

    subroutine integer_alloc_1d(vec, m, nm, errco, emsg)
        integer(8), dimension(:), allocatable, intent(inout) :: vec
        integer(8) m
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m), STAT=alloc_stat)
        call check_alloc(alloc_stat, m, nm, NBERROR_130, errco, emsg)
        end subroutine


    subroutine double_alloc_1d(vec, m, nm, errco, emsg)
        double precision, dimension(:), allocatable, intent(inout) :: vec
        integer(8) m
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m), STAT=alloc_stat)
        call check_alloc(alloc_stat, m, nm, NBERROR_130, errco, emsg)
        end subroutine

    subroutine complex_alloc_1d(vec, m, nm, errco, emsg)
        complex(8), dimension(:), allocatable, intent(inout) :: vec
        integer(8) m
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m), STAT=alloc_stat)
        call check_alloc(alloc_stat, m, nm, NBERROR_130, errco, emsg)
        end subroutine

    subroutine integer_alloc_2d(vec, m, n, nm, errco, emsg)
        integer(8), dimension(:,:), allocatable, intent(inout) :: vec
        integer(8) m, n
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m,n), STAT=alloc_stat)
        call check_alloc(alloc_stat, m*n, nm, NBERROR_130, errco, emsg)
        end subroutine


    subroutine double_alloc_2d(vec, m, n, nm, errco, emsg)
        double precision, dimension(:,:), allocatable, intent(inout) :: vec
        integer(8) m,n
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m,n), STAT=alloc_stat)
        call check_alloc(alloc_stat, m*n, nm, NBERROR_130, errco, emsg)
        end subroutine

    subroutine complex_alloc_2d(vec, m, n, nm, errco, emsg)
        complex(8), dimension(:,:), allocatable, intent(inout) :: vec
        integer(8) m,n
        character(*) nm
        integer, intent(out) :: errco
        character(len=EMSG_LENGTH), intent(out) :: emsg

        integer alloc_stat

        allocate(vec(m,n), STAT=alloc_stat)
        call check_alloc(alloc_stat, m*n, nm, NBERROR_130, errco, emsg)
        end subroutine


    end module

