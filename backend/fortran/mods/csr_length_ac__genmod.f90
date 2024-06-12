        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:47:54 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CSR_LENGTH_AC__genmod
          INTERFACE 
            SUBROUTINE CSR_LENGTH_AC(NEL,N_DDL,NEQ,NNODES,TABLE_N_E_F,  &
     &INEQ,COL_IND,ROW_PTR,NONZ_MAX,NONZ,MAX_ROW_LEN,IPOINTER,INT_MAX,  &
     &DEBUG)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: NEQ
              INTEGER(KIND=8) :: N_DDL
              INTEGER(KIND=8) :: NEL
              INTEGER(KIND=8) :: TABLE_N_E_F(NNODES,NEL)
              INTEGER(KIND=8) :: INEQ(3,N_DDL)
              INTEGER(KIND=8) :: COL_IND(*)
              INTEGER(KIND=8) :: ROW_PTR(NEQ+1)
              INTEGER(KIND=8) :: NONZ_MAX
              INTEGER(KIND=8) :: NONZ
              INTEGER(KIND=8) :: MAX_ROW_LEN
              INTEGER(KIND=8) :: IPOINTER
              INTEGER(KIND=8) :: INT_MAX
              INTEGER(KIND=8) :: DEBUG
            END SUBROUTINE CSR_LENGTH_AC
          END INTERFACE 
        END MODULE CSR_LENGTH_AC__genmod
