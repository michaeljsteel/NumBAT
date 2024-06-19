
subroutine check_alloc(alloc_stat, reqsz, nm, ec_on_fail, errco, emsg)

   use numbatmod

   character(len=*), intent(in) :: nm
   integer, intent(out) :: errco

   integer(8), intent(in) :: reqsz
   integer, intent(in)   :: ec_on_fail
   integer, intent(in) :: alloc_stat

   character(len=EMSG_LENGTH), intent(out) :: emsg

   errco = 0
   if (alloc_stat /= 0) then
      write (emsg,*) "Memory allocation for array ", nm, " failed.", &
         "  Size requested = ", reqsz , " bytes."
      errco = ec_on_fail
   endif

end
