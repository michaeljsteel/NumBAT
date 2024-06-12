        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:02 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAT_EL_ENERGY__genmod
          INTERFACE 
            SUBROUTINE MAT_EL_ENERGY(XEL,BETA,C_TENSOR_EL,MAT_P)
              REAL(KIND=8) :: XEL(2,6)
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: C_TENSOR_EL(6,6)
              COMPLEX(KIND=8) :: MAT_P(18,18)
            END SUBROUTINE MAT_EL_ENERGY
          END INTERFACE 
        END MODULE MAT_EL_ENERGY__genmod
