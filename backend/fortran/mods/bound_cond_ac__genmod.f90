        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:53 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BOUND_COND_AC__genmod
          INTERFACE 
            SUBROUTINE BOUND_COND_AC(I_COND,NPT,NEQ,TYPE_NOD,INEQ,DEBUG)
              INTEGER(KIND=8) :: NPT
              INTEGER(KIND=8) :: I_COND
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: TYPE_NOD(NPT)
              INTEGER(KIND=8) :: INEQ(3,NPT)
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE BOUND_COND_AC
          END INTERFACE 
        END MODULE BOUND_COND_AC__genmod
