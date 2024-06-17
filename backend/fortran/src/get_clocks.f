
      subroutine get_clocks(systime, cputime)
C     Returns system (wall time) in seconds, and cpu time in seconds
C     nanosec may be microsec on some systems

      integer(8) isystime
      double precision systime, cputime
      double precision nanosec

      parameter (nanosec=1.d-9)

      call system_clock(isystime)
      call cpu_time(cputime)

      systime = nanosec*isystime

      end subroutine get_clocks

