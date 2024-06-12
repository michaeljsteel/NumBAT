        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:08 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NORMALISATION__genmod
          INTERFACE 
            SUBROUTINE NORMALISATION(NVAL,NEL,NNODES,SOLN_K1,SOLN_K2,   &
     &MAT_OVERLAP)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              COMPLEX(KIND=8) :: SOLN_K1(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: SOLN_K2(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: MAT_OVERLAP(NVAL,NVAL)
            END SUBROUTINE NORMALISATION
          END INTERFACE 
        END MODULE NORMALISATION__genmod
