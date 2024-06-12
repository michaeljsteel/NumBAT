        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 12 18:53:07 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHI2_GRAD__genmod
          INTERFACE 
            SUBROUTINE PHI2_GRAD(INODE,NNODES,MAT_JAC,VEC_GRAD)
              INTEGER(KIND=8) :: NNODES
              INTEGER(KIND=8) :: INODE
              REAL(KIND=8) :: MAT_JAC(2,2)
              REAL(KIND=8) :: VEC_GRAD(2,NNODES)
            END SUBROUTINE PHI2_GRAD
          END INTERFACE 
        END MODULE PHI2_GRAD__genmod
