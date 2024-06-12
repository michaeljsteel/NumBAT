        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:51 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ARRAY_SOL_AC__genmod
          INTERFACE 
            SUBROUTINE ARRAY_SOL_AC(NUM_MODES,N_MSH_EL,N_MSH_PTS,NEQ,   &
     &NNODES,IINDEX,TABLE_NOD,TYPE_EL,INEQ,X,V_CMPLX,V_TMP,MODE_POL,    &
     &SOL_0,SOL)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: N_MSH_PTS
              INTEGER(KIND=8) :: N_MSH_EL
              INTEGER(KIND=8) :: NUM_MODES
              INTEGER(KIND=8) :: IINDEX(*)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,N_MSH_EL)
              INTEGER(KIND=8) :: TYPE_EL(N_MSH_EL)
              INTEGER(KIND=8) :: INEQ(3,N_MSH_PTS)
              REAL(KIND=8) :: X(2,N_MSH_PTS)
              COMPLEX(KIND=8) :: V_CMPLX(NUM_MODES)
              COMPLEX(KIND=8) :: V_TMP(NUM_MODES)
              COMPLEX(KIND=8) :: MODE_POL(4,NUM_MODES)
              COMPLEX(KIND=8) :: SOL_0(NEQ,NUM_MODES)
              COMPLEX(KIND=8) :: SOL(3,NNODES,NUM_MODES,N_MSH_EL)
            END SUBROUTINE ARRAY_SOL_AC
          END INTERFACE 
        END MODULE ARRAY_SOL_AC__genmod
