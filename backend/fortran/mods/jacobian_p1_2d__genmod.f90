        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:00 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE JACOBIAN_P1_2D__genmod
          INTERFACE 
            SUBROUTINE JACOBIAN_P1_2D(X,XEL,NNODES,X_G,DET_JACOBIAN,    &
     &MAT_B_0,MAT_T)
              INTEGER(KIND=8) :: NNODES
              REAL(KIND=8) :: X(2)
              REAL(KIND=8) :: XEL(2,NNODES)
              REAL(KIND=8) :: X_G(2)
              REAL(KIND=8) :: DET_JACOBIAN
              REAL(KIND=8) :: MAT_B_0(2,2)
              REAL(KIND=8) :: MAT_T(2,2)
            END SUBROUTINE JACOBIAN_P1_2D
          END INTERFACE 
        END MODULE JACOBIAN_P1_2D__genmod
