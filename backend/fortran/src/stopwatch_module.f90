module class_stopwatch
   use numbatmod

   implicit none
   private



   type, public :: stopwatch

      double precision :: cpu_t0, cpu_t1, sys_t0, sys_t1
   contains
      procedure :: reset => stopwatch_reset
      procedure :: stop => stopwatch_stop
      procedure :: report => stopwatch_report
      procedure :: cpu_time => stopwatch_cpu_time
      procedure :: sys_time => stopwatch_sys_time
      procedure :: to_string => stopwatch_to_string

   end type stopwatch

contains

   subroutine stopwatch_reset(this)
      class(stopwatch), intent(in) :: this
      call get_clocks(this%sys_t0, this%cpu_t0)
   end subroutine stopwatch_reset

   subroutine stopwatch_timecheck(this)
    class(stopwatch), intent(in) :: this
    call get_clocks(this%sys_t1, this%cpu_t1)
 end subroutine stopwatch_timecheck

   subroutine stopwatch_stop(this)
      class(stopwatch), intent(in) :: this
      call get_clocks(this%sys_t1, this%cpu_t1)
   end subroutine stopwatch_stop

   function stopwatch_cpu_time(this) result(dt_cpu)
      class(stopwatch), intent(in) :: this
      double precision :: dt_cpu
      call stopwatch_timecheck(this)
      dt_cpu = this%cpu_t1 - this%cpu_t0
   end function stopwatch_cpu_time

   function stopwatch_sys_time(this) result(dt_sys)
      class(stopwatch), intent(in) :: this
      double precision :: dt_sys
      call stopwatch_timecheck(this)
      dt_sys = this%sys_t1 - this%sys_t0
   end function stopwatch_sys_time


   function stopwatch_to_string(this) result(tstr)
      class(stopwatch), intent(in) :: this

      character(len=:), allocatable :: tstr
      integer(8),  parameter :: buflen = 512
      character(len=buflen) :: buffer

      write(buffer, '(A, f0.2, A, f0.2, A)') 'cpu time = ', this%cpu_time(), &
         ' secs, wall time = ', this%sys_time(), ' secs'
      tstr = trim(buffer)

   end function


   subroutine stopwatch_report(this, ui)
      class(stopwatch), intent(in) :: this
      integer(8), intent(in) :: ui

      write(ui, *) this%to_string()

   end subroutine stopwatch_report

end module class_stopwatch
