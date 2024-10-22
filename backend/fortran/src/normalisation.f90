!  TODO: this seems to normalise by the wrong factor: 1/|z|^2 instead of 1/|z|
!  Is this ever actually used?

subroutine normalise_fields (n_modes, n_elts, nnodes, m_evecs, mat_overlap)
   use numbatmod

   integer(8) n_modes, n_elts, nnodes
   complex(8) m_evecs(3,nnodes+7,n_modes,n_elts)
   complex(8) mat_overlap(n_modes,n_modes)

   integer(8) ival
   complex(8) z_tmp1, z_tmp2

   double precision, parameter :: min_abs = 1.0d-8

   do ival=1,n_modes
      z_tmp1 = sqrt(mat_overlap(ival,ival))

      !  Don't bother rescaling fields which are super tiny and therefore probably broken anyway
      if (abs(z_tmp1) .gt. min_abs) then
         z_tmp2 =  1.0d0/z_tmp1**2
         m_evecs(:,:,ival,:) = m_evecs(:,:,ival,:) * z_tmp2
      endif
   enddo

end subroutine

