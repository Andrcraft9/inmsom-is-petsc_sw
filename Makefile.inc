MAKE = make

# PETSc flags!
FFLAGS = -cpp -w -ffree-line-length-0 -O3 -Wtabs -fopenmp ${PETSC_FC_INCLUDES} -I./Inc
CFLAGS = ${PETSC_CC_INCLUDES}

RM = rm -rf

# PETSC options
#
RUN_PETSC_OPTS = -draw_pause -1

# -ksp_view
# -ksp_plot_eigenvalues_explicitly
PETSC_KSP_OPTS = -ksp_type gmres

#
PETSC_PC_OPTS =  -pc_type asm
