
!subroutine check_alloc(alloc_stat, reqsz, nm, ec_on_fail, errco, emsg)
subroutine check_alloc(alloc_stat, reqsz, nm, ec_on_fail, nberr)

   use numbatmod

   integer(8),  intent(in) :: alloc_stat
   integer(8), intent(in) :: reqsz
   character(len=*), intent(in) :: nm
   integer(8),  intent(in)   :: ec_on_fail
   type(NBError) :: nberr

   !integer(8),  intent(out) :: errco

   character(len=EMSG_LENGTH) :: emsg

   !errco = 0
   if (alloc_stat /= 0) then
      write (emsg,*) "Memory allocation for array ", nm, " failed.", &
         "  Size requested = ", reqsz , " bytes."
      !errco = ec_on_fail
      call nberr%set(ec_on_fail, emsg)
   endif

end
