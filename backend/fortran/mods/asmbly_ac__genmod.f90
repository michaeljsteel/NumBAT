        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:52 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASMBLY_AC__genmod
          INTERFACE 
            SUBROUTINE ASMBLY_AC(I_BASE,NEL,NPT,NEQ,NNODES,SHIFT,BETA,  &
     &NB_TYP_EL,RHO,C_TENSOR,TABLE_NOD,TYPE_EL,INEQ,X,NONZ,ROW_IND,     &
     &COL_PTR,MAT1_RE,MAT1_IM,MAT2,I_WORK,SYMMETRY_FLAG,DEBUG)
              INTEGER(KIND=8) :: NONZ
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: I_BASE
              COMPLEX(KIND=8) :: SHIFT
              COMPLEX(KIND=8) :: BETA
              COMPLEX(KIND=8) :: RHO(NB_TYP_EL)
              COMPLEX(KIND=8) :: C_TENSOR(6,6,NB_TYP_EL)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              INTEGER(KIND=8) :: INEQ(3,NPT)
              REAL(KIND=8) :: X(2,NPT)
              INTEGER(KIND=8) :: ROW_IND(NONZ)
              INTEGER(KIND=8) :: COL_PTR(NEQ+1)
              REAL(KIND=8) :: MAT1_RE(NONZ)
              REAL(KIND=8) :: MAT1_IM(NONZ)
              COMPLEX(KIND=8) :: MAT2(NONZ)
              INTEGER(KIND=8) :: I_WORK(3*NPT)
              INTEGER(KIND=8) :: SYMMETRY_FLAG
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE ASMBLY_AC
          END INTERFACE 
        END MODULE ASMBLY_AC__genmod
