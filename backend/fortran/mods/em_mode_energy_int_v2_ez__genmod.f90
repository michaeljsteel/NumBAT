        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:57 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EM_MODE_ENERGY_INT_V2_EZ__genmod
          INTERFACE 
            SUBROUTINE EM_MODE_ENERGY_INT_V2_EZ(K_0,NVAL,NEL,NPT,       &
     &NNODES_P2,TABLE_NOD,X,BETAS,SOLN_K1,OVERLAP)
              INTEGER(KIND=8) :: NNODES_P2
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              REAL(KIND=8) :: K_0
              INTEGER(KIND=8) :: TABLE_NOD(NNODES_P2,NEL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: BETAS(NVAL)
              COMPLEX(KIND=8) :: SOLN_K1(3,NNODES_P2+7,NVAL,NEL)
              COMPLEX(KIND=8) :: OVERLAP(NVAL)
            END SUBROUTINE EM_MODE_ENERGY_INT_V2_EZ
          END INTERFACE 
        END MODULE EM_MODE_ENERGY_INT_V2_EZ__genmod
