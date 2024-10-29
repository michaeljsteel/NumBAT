
! Construction of the matrix of power flow for one elt
!	(integral of the z-component of the acoustic Poynting vector)

! p[i,j] = \int_{one elt} \dx\dy  c_zjkl u_j^* d_k u_l
! (Basis functions are real so star on u_j^* is ignored)
! Here we perform the sum over k, but not for j and l which
!   run over basis funcs: [1x, 1y, 1z, 2x, 2y, 2z, ..., 6x, 6y, 6z]
!
! So  p[i,j] = \int_{one elt} \dx\dy (
!                 c_zjxl u_j (d_x u_l) + c_zjyl u_j d_y u_l + c_zjyl u_j (i beta) u_l

! Only include terms which will give an imaginary contribution
! (since there is a factor (2 i Omega) in the calling routine)
! This means the cases are
!    -  u_x/y and u_x/y with a z derivative
!    -  u_x/y and u_z with an x/y derivative
!    -  u_z and u_z with a z derivative
!   making seven cases
!
! Below we write u for u_j, and v for u_l

! should be named mat_el_powerflow
subroutine mat_el_powerflow (q_AC, stiff_C_IJ, basfuncs, mat_P)

    use numbatmod
    use class_BasisFunctions

    complex(8) q_AC
    complex(8) mat_P(18,18)
    complex(8) stiff_C_IJ(6,6)
    type(BasisFunctions) basfuncs


    ! Locals

    double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)

    double precision mat_T_tr(2,2)
    double precision det_b

    complex(8) z_tmp1
    integer(8) i, j, i_p, j_p, i_xyz,  j_xyz
    integer(8) debug

    !    Compute the Affine mappings from the current triangle to the
    !     reference unit triangle.
    !    Integration will be performed on the reference unit triangle

    det_b = basfuncs%det

    !	mat_T_tr = Tanspose(mat_T)

    mat_T_tr(1,1) = basfuncs%mat_T(1,1)
    mat_T_tr(1,2) = basfuncs%mat_T(2,1)
    mat_T_tr(2,1) = basfuncs%mat_T(1,2)
    mat_T_tr(2,2) = basfuncs%mat_T(2,2)

    call find_overlaps_p2_p2(p2_p2, det_b)
    call find_overlaps_p2_p2x (p2_p2x, mat_T_tr, det_b)
    call find_overlaps_p2_p2y (p2_p2y, mat_T_tr, det_b)


    mat_P = D_ZERO

    !=================  Construction of the matrix of power flow   ==================
    !                   (integral of the z-component of the acoustic Poynting vector)

    do i=1,P2_NODES_PER_EL  ! Loop over basis funcs
       do i_xyz=1,3
          i_p = 3*(i-1) + i_xyz

          do j=1,P2_NODES_PER_EL  ! Loop over basis funcs
             do j_xyz=1,3
                j_p = 3*(j-1) + j_xyz

                if (i_xyz == 1 .and. j_xyz == 1) then
                   ! Agree
                   ! c_zxkx ux d_k vx -> c_zxzx ux i q vx = i q C(5,5)  ux vx
                   z_tmp1 = C_IM_ONE  * q_AC *  stiff_C_IJ(5,5) * p2_p2(i,j)

                elseif (i_xyz == 1 .and. j_xyz == 3) then
                   ! DISAGREE
                   ! c_zxkz ux d_k vz ->
                   !  = c_zxxx ux vz_x + c_zxyx ux vz_y
                   !  = c_51 ux vz_x + c_56 ux vz_y
                   ! z_tmp1 = stiff_C_IJ(5,1) * p2_p2x(i,j) + stiff_C_IJ(5,6) * p2_p2y(i,j)

                   ! Orig:  C(5,5) * u_x * S_xz
                   z_tmp1 = p2_p2x(i,j) * stiff_C_IJ(5,5)

                elseif (i_xyz == 2 .and. j_xyz == 2) then
                   ! Agree
                   ! c_zyky uy d_k vy ->  c_zyzy uy vy_z  = i q c_44 uy vy

                   !   C(4,4) * u_x * S_zy
                   z_tmp1 = C_IM_ONE * q_AC * stiff_C_IJ(4,4) * p2_p2(i,j)

                elseif (i_xyz == 2 .and. j_xyz == 3) then
                   ! DISAGREE
                   ! c_zykz uy d_k vz -> c_zyxz uy vz_x + c_zyyz uy vz_y
                   !  = c_45 uy vz_x + c_44 uy vz_y
                   ! z_tmp1 = stiff_C_IJ(4,5) * p2_p2x(i,j) + stiff_C_IJ(4,4) * p2_p2y(i,j)

                   !  Orig: C(4,4) * u_x * S_yz
                   z_tmp1 = p2_p2y(i,j) * stiff_C_IJ(4,4)

                elseif (i_xyz == 3 .and. j_xyz == 1) then
                   ! DISAGREE
                   ! c_zzkx uz d_k vx ->
                   !  = c_zzxx uz vx_x + c_zzyx uz vx_y
                   !  = c_31 uz vx_x + c_36 uz vx_y
                   ! z_tmp1 = stiff_C_IJ(3,1) * p2_p2x(i,j) + stiff_C_IJ(3,6) * p2_p2y(i,j)

                   !  Orig: C(3,1) * u_x * S_xx
                   z_tmp1 =  stiff_C_IJ(3,1) * p2_p2x(i,j)

                elseif (i_xyz == 3 .and. j_xyz == 2) then
                   ! DISAGREE
                   ! c_zzky uz d_k vy ->
                   !  = c_zzxy uz vy_x + c_zzyy uz vy_y
                   !  = c_36 uz vy_x + c_32 uz vy_y
                   ! z_tmp1 = stiff_C_IJ(3,6) * p2_p2x(i,j) + stiff_C_IJ(3,2) * p2_p2y(i,j)

                   !  Orig: C(3,2) * u_x * S_yy
                   z_tmp1 = p2_p2y(i,j) * stiff_C_IJ(3,2)

                elseif (i_xyz == 3 .and. j_xyz == 3) then
                   ! Agree
                   ! c_zzkz uz d_k vz ->  c_zzzz iq uz vz_z  = i q c_33 uz vz

                   z_tmp1 = C_IM_ONE * q_AC *  stiff_C_IJ(3,3) * p2_p2(i,j)

                else

                   z_tmp1 = D_ZERO

                endif

                mat_P(i_p,j_p) = mat_P(i_p,j_p)  + z_tmp1

             enddo
          enddo
       enddo
    enddo

 end subroutine mat_el_powerflow
