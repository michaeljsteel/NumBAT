module class_Stopwatch
    use numbatmod

    implicit none
    private



    type, public :: Stopwatch

    double precision :: cpu_t0, cpu_t1, sys_t0, sys_t1
    contains
    procedure :: reset => stopwatch_reset
    procedure :: stop => stopwatch_stop
    procedure :: report => stopwatch_report
    procedure :: cpu_time => stopwatch_cpu_time
    procedure :: sys_time => stopwatch_sys_time
    procedure :: to_string => stopwatch_to_string

    end type Stopwatch

    contains

    subroutine stopwatch_reset(this)
        class(Stopwatch), intent(in) :: this
        call get_clocks(this%sys_t0, this%cpu_t0)
        end subroutine stopwatch_reset

        subroutine stopwatch_stop(this)
            class(Stopwatch), intent(in) :: this
            call get_clocks(this%sys_t1, this%cpu_t1)
            end subroutine stopwatch_stop

            function stopwatch_cpu_time(this) result(dt_cpu)
                class(Stopwatch), intent(in) :: this
                double precision :: dt_cpu
                dt_cpu = this%cpu_t1 - this%cpu_t0
            end function stopwatch_cpu_time

            function stopwatch_sys_time(this) result(dt_sys)
                class(Stopwatch), intent(in) :: this
                double precision :: dt_sys
                dt_sys = this%sys_t1 - this%sys_t0
            end function stopwatch_sys_time


            function stopwatch_to_string(this) result(tstr)
                class(Stopwatch), intent(in) :: this

            character(len=:), allocatable :: tstr
            integer, parameter :: buflen = 512
            character(len=buflen) :: buffer

            write(buffer, '(A, f0.2, A, f0.2, A)') 'cpu time = ', this%cpu_time(), &
                    ' secs, wall time = ', this%sys_time(), ' secs'
            tstr = trim(buffer)

            end function


            subroutine stopwatch_report(this, ui)
                class(Stopwatch), intent(in) :: this
                integer(8), intent(in) :: ui

                call this%stop()
                write(ui, *) this%to_string()

                end subroutine stopwatch_report

    end module class_Stopwatch