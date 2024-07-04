
subroutine conv_gmsh(geo_name, assertions_on, errco, emsg)

   use numbatmod

   !TODO: f2py doesn't seem to be able to handle this length being set by the #include and preprocessor
   integer :: assertions_on
   character(len=*), intent(in) :: geo_name

   integer, intent(out) :: errco
   character(len=EMSG_LENGTH), intent(out) :: emsg

   errco = 0
   call conv_gmsh_impl(geo_name, assertions_on, errco, emsg)

end
