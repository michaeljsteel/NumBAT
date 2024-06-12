        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:55 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EM_MODE_E_ENERGY_INT__genmod
          INTERFACE 
            SUBROUTINE EM_MODE_E_ENERGY_INT(NVAL,NEL,NPT,NNODES,        &
     &TABLE_NOD,TYPE_EL,NB_TYP_EL,N_LST,X,SOLN_EM,OVERLAP)
              INTEGER(KIND=8) :: NB_TYP_EL
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: NVAL
              INTEGER(KIND=8) :: TABLE_NOD(NNODES,NEL)
              INTEGER(KIND=8) :: TYPE_EL(NEL)
              COMPLEX(KIND=8) :: N_LST(NB_TYP_EL)
              REAL(KIND=8) :: X(2,NPT)
              COMPLEX(KIND=8) :: SOLN_EM(3,NNODES+7,NVAL,NEL)
              COMPLEX(KIND=8) :: OVERLAP(NVAL)
            END SUBROUTINE EM_MODE_E_ENERGY_INT
          END INTERFACE 
        END MODULE EM_MODE_E_ENERGY_INT__genmod
