include Makefile.inc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SRCCONTROL =	\
	Control/dyncall.f90      \
	Control/forc_atm.f90     \
	Control/init_arrays.f90  \
	Control/init_pars.f90    \
	Control/output.f90

SRCSERVICE =	    \
	Service/atm2oc_interpol.f90       \
	Service/basin_parameters.f90      \
	Service/bc_time_interpol.f90      \
	Service/grid_construction.f90     \
	Service/input_output_data.f90     \
	Service/oc_algorithms.f90         \
	Service/read_write_parameters.f90 \
	Service/rw_ctl_file.f90           \
	Service/time_tools.f90


SRCFUNCTION =	\
	Function/depth.f90   \
	Function/mixing.f90  \
	Function/vel_ssh.f90 \
	Function/fluxes.f90


SRCMODULES = 	\
	Modules/mod_atm2oc_interpol.f90      \
	Modules/mod_atm_pars.f90             \
	Modules/mod_atm_variables.f90        \
	Modules/mod_basin_grid.f90           \
	Modules/mod_main_basin_pars.f90      \
	Modules/mod_mpi_parallel_tools.f90   \
	Modules/mod_ocean_variables.f90      \
	Modules/mod_rec_length.f90           \
	Modules/mod_time_integration.f90     \
	Modules/parallel_ssh_solve.f90  

LIBS = 	\
	$(LIBAMD)/amd.f


#inmsom:
#order is important
#	$(FC) -o inmsom $(FCFLAGS) $(SRCMODULES) inmsom_head.f90  $(SRCFUNCTION) $(SRCCONTROL) $(SRCSERVICE)  \
#	-L $(LIBMAINDIR)/$(LIBSPARSKIT) -lskit -L $(LIBMAINDIR)/$(LIBBLAS) -lblas $(LIBS)

inmsom: 
	-${FLINKER} -o inmsom $(SRCMODULES) inmsom_head.f90  $(SRCFUNCTION) $(SRCCONTROL) $(SRCSERVICE) ${PETSC_KSP_LIB} 

clean_:
	$(RM) *.o *.mod

run_seq:
	-@${MPIEXEC} -n 1 ./inmsom -ksp_type gmres -pc_type none

#-ksp_type gmres
run:
	-@${MPIEXEC} -n 4 ./inmsom -ksp_type gmres -pc_type none -ksp_view -log_summary


