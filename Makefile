########################################################################
#
#                           MAKEFILE FOR ADC                         
#
########################################################################

#-----------------------------------------------------------------------
# Inclusion of other makefiles
#-----------------------------------------------------------------------
include ${SLEPC_DIR}/conf/slepc_common

#-----------------------------------------------------------------------
# Compiler flags
#
# N.B. we now have to use the C preprocessor due to interfacing with
# PETSc/SLEPc
#-----------------------------------------------------------------------
F90	= gfortran
F77	= gfortran
CC	= gcc

F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O2
CCOPTS  = -g -O0

#-----------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------
LIBS= -L/usr/lib64 ${LIB_LAPACK} ${LIB_BLAS}

SLEPC_LIBS=${SLEPC_EPS_LIB} -I${PETSC_DIR}/include -I${SLEPC_DIR}/include

#-----------------------------------------------------------------------
# Define object files
#-----------------------------------------------------------------------
# Main ADC code
MULTI = accuracy.o \
	printing.o \
	timer.o \
	lapack.o \
	dgefa.o \
	dgedi.o \
	math.o \
	matrix_tools.o \
	os_integral_operators.o \
	gamess_internal.o \
	import_gamess.o \
	integral_tools.o \
	integrals_mo2e.o \
        diis.o \
        sort_tools.o \
        block_diag.o \
        biorthogonal_tools.o \
        scf_tools.o \
        fock_tools.o 

ADC =   constants.o \
	parameters.o \
	misc.o \
	load_electronic_structure.o \
        scf_electronic_structure.o \
	external_diag.o \
	filetools.o \
	adc_ph.o \
	dipole_ph.o \
	D_matrix.o \
	read_param.o \
	sym_allowed_exc.o \
	select_fano.o \
	get_matrix.o \
	get_matrix_dipole_complete.o \
	get_moment.o fspacetrial.o \
	fspace2.o \
	davmod.o \
	band_lanczos.o \
	photoionisation.o \
	Propagate.o \
        guessvecs.o \
	master_adc1_prop.o \
	master_adc2_prop.o \
	master_adc2ext_prop.o \
	master_adc2_ener.o \
	main_draft1.o

OBJECTS = $(MULTI) $(ADC)

# Stieltjes imaging code
STIELTJES = stieltjes/qmath.o \
	stieltjes/pythag_quad.o \
	stieltjes/tql2_quad.o \
        stieltjes/globalmod.o \
	stieltjes/stieltjes_s3_modified.o \
	stieltjes/main_stieltjes.o 

STIELTJES_OBJ = qmath.o \
	pythag_quad.o \
	tql2_quad.o \
        globalmod.o \
	stieltjes_s3_modified.o \
	main_stieltjes.o 

#-----------------------------------------------------------------------
# Rules to create the program
#-----------------------------------------------------------------------
adc: $(OBJECTS)
	$(F90) $(F90OPTS) $(OBJECTS) $(LIBS) $(SLEPC_LIBS) -o adc.x 

stieltjes: $(STIELTJES)
	$(F90) $(F90OPTS) $(STIELTJES_OBJ) -o stieltjes.x

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
