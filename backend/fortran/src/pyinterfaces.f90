
!#include "numbat_decl.h"


subroutine conv_gmsh(geoname, assertions_on, errco, emsg)

   use numbatmod

   !TODO: f2py doesn't seem to be able to handle this length being set by the #include and preprocessor
   integer :: assertions_on, errco

   character(len=FNAME_LENGTH) :: geoname
   character(len=EMSG_LENGTH) :: emsg

!f2py intent(in) geoname, assertions_on
!f2py intent(out) errco, emsg

   errco = 0
   call conv_gmsh_impl(geoname, assertions_on, errco, emsg)

end
