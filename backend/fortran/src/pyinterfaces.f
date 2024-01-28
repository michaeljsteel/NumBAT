C  Interface file for f2py
C  Intended to reduce the number of symbols exposed to python within NumBAT.so      
   
!#include "numbat_decl.h"


      subroutine conv_gmsh(geoname, assertions_on, errco, emsg)

      ! use numbatmod
      implicit none
 
      !TODO: f2py doesn't seem to be able to handle this length being set by the #include and preprocessor 
      integer, parameter :: STRINGLEN =  1024   
      character geoname*(STRINGLEN)
      integer assertions_on, errco
      character emsg*(STRINGLEN)

Cf2py intent(in) geoname, assertions_on
Cf2py intent(out) errco, emsg

      errco = 0
      call conv_gmsh_impl(geoname, assertions_on, errco, emsg)
      
      end