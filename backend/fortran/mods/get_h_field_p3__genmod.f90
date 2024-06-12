        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:58 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_H_FIELD_P3__genmod
          INTERFACE 
            SUBROUTINE GET_H_FIELD_P3(NNODES_P2,K_0,BETA1,MAT_T,        &
     &E_FIELD_EL,EZ_FIELD_EL_P3,H_FIELD_EL)
              INTEGER(KIND=8) :: NNODES_P2
              REAL(KIND=8) :: K_0
              COMPLEX(KIND=8) :: BETA1
              REAL(KIND=8) :: MAT_T(2,2)
              COMPLEX(KIND=8) :: E_FIELD_EL(3,NNODES_P2)
              COMPLEX(KIND=8) :: EZ_FIELD_EL_P3(10)
              COMPLEX(KIND=8) :: H_FIELD_EL(3,NNODES_P2)
            END SUBROUTINE GET_H_FIELD_P3
          END INTERFACE 
        END MODULE GET_H_FIELD_P3__genmod
