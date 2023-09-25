      subroutine check_alloc(stat, reqsz, nm, ec, errco, emsg)

      implicit none

      character*2048 nm, emsg
      integer*8      ec, errco
      integer :: stat, reqsz

      errco = 0
      if (stat .ne. 0) then
          write (emsg,*) "Memory allocation for ", nm, "failed.",
     *          "Size requested = ", reqsz , " bytes."
          errco = ec
      endif

      end
