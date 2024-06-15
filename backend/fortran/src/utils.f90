
subroutine check_alloc(stat, reqsz, nm, ec, errco, emsg)

   use numbatmod

   character(len=*), intent(in) :: nm


   integer*8 :: reqsz
   integer   :: ec, errco 
   integer :: stat

   character(len=EMSG_LENGTH), intent(out) :: emsg

   errco = 0
   if (stat .ne. 0) then
      write (emsg,*) "Memory allocation for ", nm, "failed.", &
         "Size requested = ", reqsz , " bytes."
      errco = ec
   endif

end
