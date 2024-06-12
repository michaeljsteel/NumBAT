        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BASIS_VEC__genmod
          INTERFACE 
            SUBROUTINE BASIS_VEC(I_EQ,I_DDL,BASIS_LIST,P2_LIST,         &
     &GRAD_P1_MAT,GRAD_P2_MAT,VEC_PHI,CURL_T_PHI)
              INTEGER(KIND=8) :: I_EQ
              INTEGER(KIND=8) :: I_DDL
              INTEGER(KIND=8) :: BASIS_LIST(4,3,4)
              REAL(KIND=8) :: P2_LIST(6)
              REAL(KIND=8) :: GRAD_P1_MAT(2,3)
              REAL(KIND=8) :: GRAD_P2_MAT(2,6)
              REAL(KIND=8) :: VEC_PHI(2)
              REAL(KIND=8) :: CURL_T_PHI
            END SUBROUTINE BASIS_VEC
          END INTERFACE 
        END MODULE BASIS_VEC__genmod
