MAKE = make

# PETSc flags!
FFLAGS = -cpp -w -ffree-line-length-0 -O3 -Wtabs -fopenmp ${PETSC_FC_INCLUDES} -I./Inc
CFLAGS = ${PETSC_CC_INCLUDES}

RM = rm -rf


