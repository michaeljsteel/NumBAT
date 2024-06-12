        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:03 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAT_EL_V3__genmod
          INTERFACE 
            SUBROUTINE MAT_EL_V3(XEL,BETA,C_TENSOR_EL,RHO_EL,MAT_K,MAT_M&
     &)
              REAL(KIND=8) :: XEL(2,6)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C_TENSOR_EL(6,6)
              COMPLEX(KIND=8) :: RHO_EL
              COMPLEX(KIND=8) :: MAT_K(18,18)
              COMPLEX(KIND=8) :: MAT_M(18,18)
            END SUBROUTINE MAT_EL_V3
          END INTERFACE 
        END MODULE MAT_EL_V3__genmod
