        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:55 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EM_MODE_ENERGY_INT__genmod
          INTERFACE 
            SUBROUTINE EM_MODE_ENERGY_INT(K_0,NVAL,NEL,NPT,NNODES,      &
     &TABLE_NOD,X,BETAS,SOLN_K1,OVERLAP)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              REAL(KIND=8) :: K_0
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: BETAS(NVAL)
              COMPLEX(KIND=8) :: SOLN_K1(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: OVERLAP(NVAL)
            END SUBROUTINE EM_MODE_ENERGY_INT
          END INTERFACE 
        END MODULE EM_MODE_ENERGY_INT__genmod
