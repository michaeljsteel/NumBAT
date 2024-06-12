        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ORTHOGONAL__genmod
          INTERFACE 
            SUBROUTINE ORTHOGONAL(NVAL,NEL,NPT,NNODES,NB_TYP_EL,PP,     &
     &TABLE_NOD,TYPE_EL,X,BETA1,BETA2,SOLN_K1,SOLN_K2,MAT_OVERLAP,      &
     &OVERLAP_FILE,PRINTALL,PAIR_WARNING,K_0)
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              COMPLEX(KIND=8) :: PP(NB_TYP_EL)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: BETA1(NVAL)
              COMPLEX(KIND=8) :: BETA2(NVAL)
              COMPLEX(KIND=8) :: SOLN_K1(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: SOLN_K2(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: MAT_OVERLAP(NVAL,NVAL)
              CHARACTER(LEN=100) :: OVERLAP_FILE
              INTEGER(KIND=8) :: PRINTALL
              INTEGER(KIND=8) :: PAIR_WARNING
              REAL(KIND=8) :: K_0
            END SUBROUTINE ORTHOGONAL
          END INTERFACE 
        END MODULE ORTHOGONAL__genmod
