#include "numbat_decl.h"

subroutine check_alloc(stat, reqsz, nm, ec, errco, emsg)

   use numbatmod

   character(len=EMSG_LENGTH) emsg

   character*2048 :: nm
   integer*8      :: ec, errco, reqsz
   integer :: stat

   errco = 0
   if (stat .ne. 0) then
      write (emsg,*) "Memory allocation for ", nm, "failed.", &
         "Size requested = ", reqsz , " bytes."
      errco = ec
   endif

end
