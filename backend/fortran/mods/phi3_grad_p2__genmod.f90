        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:08 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHI3_GRAD_P2__genmod
          INTERFACE 
            SUBROUTINE PHI3_GRAD_P2(INODE,NNODES_P3,MAT_JAC,VEC_GRAD)
              INTEGER(KIND=8) :: NNODES_P3
              INTEGER(KIND=8) :: INODE
              REAL(KIND=8) :: MAT_JAC(2,2)
              REAL(KIND=8) :: VEC_GRAD(2,NNODES_P3)
            END SUBROUTINE PHI3_GRAD_P2
          END INTERFACE 
        END MODULE PHI3_GRAD_P2__genmod
