        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:48:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MODE_ENERGY__genmod
          INTERFACE 
            SUBROUTINE MODE_ENERGY(NVAL,NEL,NPT,NNODES,N_CORE,TABLE_NOD,&
     &TYPE_EL,NB_TYP_EL,EPS_EFF,X,SOL,BETA1,MODE_POL)
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              INTEGER(KIND=8) :: N_CORE(2)
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              COMPLEX(KIND=8) :: EPS_EFF(NB_TYP_EL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: SOL(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: BETA1(NVAL)
              COMPLEX(KIND=8) :: MODE_POL(4,NVAL)
            END SUBROUTINE MODE_ENERGY
          END INTERFACE 
        END MODULE MODE_ENERGY__genmod
