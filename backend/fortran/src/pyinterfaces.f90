
subroutine conv_gmsh(geoname, assertions_on, errco, emsg)

   use numbatmod

   !TODO: f2py doesn't seem to be able to handle this length being set by the #include and preprocessor
   integer :: assertions_on
   character(len=FNAME_LENGTH) :: geoname

   integer, intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   errco = 0
   call conv_gmsh_impl(geoname, assertions_on, errco, emsg)

end
