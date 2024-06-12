        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 19:15:43 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALC_EM_MODES__genmod
          INTERFACE 
            SUBROUTINE CALC_EM_MODES(N_MODES,LAMBDA,DIMSCALE_IN_M,      &
     &BLOCH_VEC,SHIFT_KSQR,E_H_FIELD,BND_CDN_I,ITERMAX,DEBUG,MESH_FILE, &
     &N_MSH_PTS,N_MSH_EL,N_TYP_EL,V_REFINDEX_N,BETA1,SOL1,MODE_POL,     &
     &TABLE_NOD,TYPE_EL,TYPE_NOD,MESH_XY,LS_MATERIAL,ERRCO,EMSG)
              INTEGER(KIND=8) :: N_TYP_EL
              INTEGER(KIND=8) :: N_MSH_EL
              INTEGER(KIND=8) :: N_MSH_PTS
              INTEGER(KIND=8) :: N_MODES
              REAL(KIND=8) :: LAMBDA
              REAL(KIND=8) :: DIMSCALE_IN_M
              REAL(KIND=8) :: BLOCH_VEC(2)
              COMPLEX(KIND=8) :: SHIFT_KSQR
              INTEGER(KIND=8) :: E_H_FIELD
              INTEGER(KIND=8) :: BND_CDN_I
              INTEGER(KIND=8) :: ITERMAX
              INTEGER(KIND=8) :: DEBUG
              CHARACTER(LEN=1000) :: MESH_FILE
              COMPLEX(KIND=8) :: V_REFINDEX_N(N_TYP_EL)
              COMPLEX(KIND=8) ,TARGET :: BETA1(N_MODES)
              COMPLEX(KIND=8) ,TARGET :: SOL1(3,13,N_MODES,N_MSH_EL)
              COMPLEX(KIND=8) :: MODE_POL(4,N_MODES)
              INTEGER(KIND=8) :: TABLE_NOD(6,N_MSH_EL)
              INTEGER(KIND=8) :: TYPE_EL(N_MSH_EL)
              INTEGER(KIND=8) :: TYPE_NOD(N_MSH_PTS)
              REAL(KIND=8) :: MESH_XY(2,N_MSH_PTS)
              COMPLEX(KIND=8) :: LS_MATERIAL(1,13,N_MSH_EL)
              INTEGER(KIND=8) :: ERRCO
              CHARACTER(LEN=2048) :: EMSG
            END SUBROUTINE CALC_EM_MODES
          END INTERFACE 
        END MODULE CALC_EM_MODES__genmod
