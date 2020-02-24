########################################################################
#
#
#                           MAKEFILE FOR ADC                         
#
########################################################################

#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------

#
# gfortran
#
#F90	= gfortran
#F77	= gfortran
#CC	= gcc
##F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fcheck=bounds -fcheck=all -fopenmp -O3 -fbacktrace
#F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O3 -fbacktrace
#CCOPTS  = -g -O0

#
# intel
#
F90	 = ifort
F77	 = ifort
CC	 = icc
F90OPTS = -cpp -g -assume byterecl -free -fopenmp -traceback -O2 -unroll0 -diag-disable 8290 -diag-disable 8291
CCOPTS  = -g -O0

#-----------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------
LIBS= ${LIB_LAPACK} ${LIB_BLAS}

#-----------------------------------------------------------------------
# Define object files
#-----------------------------------------------------------------------

########################################################################
# ADC code
########################################################################
MULTI = source/multi/accuracy.o \
	source/multi/printing.o \
	source/multi/timer.o \
	source/multi/lapack.o \
	source/multi/dgefa.o \
	source/multi/dgedi.o \
	source/multi/math.o \
	source/multi/matrix_tools.o \
	source/multi/os_integral_operators.o \
	source/multi/gamess_internal.o \
	source/multi/import_gamess.o \
	source/multi/integral_tools.o \
	source/multi/integrals_mo2e.o \
        source/multi/diis.o \
        source/multi/sort_tools.o \
        source/multi/block_diag.o \
        source/multi/biorthogonal_tools.o \
        source/multi/scf_tools.o \
        source/multi/fock_tools.o 

INCLUDE=source/include/constants.o \
	source/include/parameters.o \
	source/include/channels.o

UTILITIES=source/utilities/timingmod.o \
	source/utilities/misc.o \
	source/utilities/lineq.o \
	source/utilities/eigchess.o \
	source/utilities/gammainc.o

IOMODULES=source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/iomodules/rdinput.o \
	source/iomodules/read_param.o \
	source/iomodules/molden.o

QCLIB=	source/qclib/vpqrsmod.o \
	source/qclib/rungamess.o \
	source/qclib/load_electronic_structure.o \
	source/qclib/scf_electronic_structure.o 

PROPAGATION= source/propagation/tdselib.o \
	source/propagation/sillib.o \
	source/propagation/csillib.o \
	source/propagation/bslib.o \
	source/propagation/rkf45rlxlib.o \
	source/propagation/specbounds.o \
	source/propagation/relaxation.o \
	source/propagation/fvecprop.o \
	source/propagation/chebyspec.o \
	source/propagation/fdstates_tdfd.o \
	source/propagation/fdstates_cfd.o \
	source/propagation/flux.o \
	source/propagation/proplib_adc2.o \
	source/propagation/proplib_adc1.o

CAP=	source/cap/auto_cap_box.o \
	source/cap/monomial_analytic.o \
	source/cap/lebedev.o \
	source/cap/atoms.o \
	source/cap/molecular_grid.o \
	source/cap/rf_cap.o \
	source/cap/basis_cap.o \
	source/cap/cap_mobas.o \
	source/cap/theta_mobas.o

ADCLIB= source/adclib/defaults.o \
	source/adclib/orbindx.o \
	source/adclib/filetools.o \
	source/adclib/adc_ph.o \
	source/adclib/dipole_ph.o \
	source/adclib/D_matrix.o \
	source/adclib/sym_allowed_exc.o \
	source/adclib/select_fano.o \
	source/adclib/get_matrix.o \
	source/adclib/get_matrix_dipole_complete.o \
	source/adclib/get_moment.o \
	source/adclib/fspacetrial.o \
	source/adclib/fspace2.o \
        source/adclib/guessvecs.o \
	source/adclib/mp2.o \
	source/adclib/electron_density.o \
	source/adclib/density_matrix.o \
	source/adclib/dyson_calc.o \
	source/adclib/dyson_io.o \
	source/adclib/target_matching.o \
	source/adclib/nto.o

EIGEN= source/eigen/block_davidson.o \
	source/eigen/dmatvec_davidson.o \
        source/eigen/block_lanczos.o \
	source/eigen/power.o

ADCCOMMON= source/adclib/adc_common.o

ADC_MAIN=source/adc/adc1_opa.o \
	source/adc/adc1_ener.o \
	source/adc/adc2_opa.o \
	source/adc/adc2_ener.o \
	source/adc/adc2_dyson.o \
	source/adc/adc2_rixs.o \
	source/adc/adc2_tpa.o \
	source/adc/adc2_autospec.o \
	source/adc/adc2_fdstates.o \
	source/adc/adc2_propagate.o \
	source/adc/adc1_propagate.o \
	source/adc/adc.o

ADC =   $(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(QCLIB) \
	$(ADCLIB) \
	$(EIGEN) \
	$(PROPAGATION) \
	$(ADCCOMMON) \
	$(CAP) \
	$(ADC_MAIN)

OBJECTS = $(MULTI) $(ADC)

ADC_OBJ=accuracy.o \
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
        fock_tools.o \
        constants.o \
	parameters.o \
	channels.o \
	misc.o \
	lineq.o \
	eigchess.o \
	gammainc.o \
	iomod.o \
	parsemod.o \
	rdinput.o \
	molden.o \
	timingmod.o \
	defaults.o \
	orbindx.o \
	vpqrsmod.o \
	rungamess.o \
	load_electronic_structure.o \
        scf_electronic_structure.o \
	filetools.o \
	adc_ph.o \
	dipole_ph.o \
	D_matrix.o \
	read_param.o \
	sym_allowed_exc.o \
	select_fano.o \
	get_matrix.o \
	get_matrix_dipole_complete.o \
	get_moment.o \
	fspacetrial.o \
	fspace2.o \
	block_davidson.o \
        dmatvec_davidson.o \
	block_lanczos.o \
	power.o \
	tdselib.o \
	sillib.o \
	csillib.o \
	bslib.o \
	rkf45rlxlib.o \
	specbounds.o \
	relaxation.o \
	chebyspec.o \
	fvecprop.o \
	fdstates_tdfd.o \
	fdstates_cfd.o \
	flux.o \
	proplib_adc2.o \
	proplib_adc1.o \
	mp2.o \
	electron_density.o \
	density_matrix.o \
	dyson_calc.o \
	dyson_io.o \
	target_matching.o \
	adc_common.o \
	nto.o \
	guessvecs.o \
	auto_cap_box.o \
	monomial_analytic.o \
	lebedev.o \
	atoms.o \
	molecular_grid.o \
	rf_cap.o \
	basis_cap.o \
	cap_mobas.o \
	theta_mobas.o \
	adc1_opa.o \
	adc1_ener.o \
	adc2_opa.o \
	adc2_ener.o \
	adc2_dyson.o \
	adc2_rixs.o \
	adc2_tpa.o \
	adc2_autospec.o \
	adc2_fdstates.o \
	adc2_propagate.o \
	adc1_propagate.o \
	adc.o

########################################################################
# Stieltjes imaging code
########################################################################
STIELTJES_AP = source/mpfun/second.o \
	source/mpfun/mpfuna.o \
	source/mpfun/mpfunbq.o \
	source/mpfun/mpfunc.o \
	source/mpfun/mpfund.o \
	source/mpfun/mpfune.o \
	source/mpfun/mpfunfq1.o \
	source/mpfun/mpmodule.o \
	source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/stieltjes_ap/pythag.o \
	source/analysis/stieltjes_ap/tql2.o \
        source/analysis/stieltjes_ap/simod.o \
        source/analysis/stieltjes_ap/stieltjes_ap.o 

STIELTJES_AP_OBJ = second.o \
	mpfuna.o \
	mpfunbq.o \
	mpfunc.o \
	mpfund.o \
	mpfune.o \
	mpfunfq1.o \
	mpmodule.o \
	constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	pythag.o \
	tql2.o \
        simod.o \
        stieltjes_ap.o

########################################################################
# Monotonicity-constrained spline interpolation code
########################################################################
MCSPLINE = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/mcspline/mcspmod.o \
	source/analysis/mcspline/mcspline.o

MCSPLINE_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	mcspmod.o \
	mcspline.o

########################################################################
# Numerical Hessian code
########################################################################
NUMHESS = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/numhess/hessmod.o \
	source/analysis/numhess/prepmod.o \
	source/analysis/numhess/calcmod.o \
	source/analysis/numhess/numhess.o

NUMHESS_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	hessmod.o \
	prepmod.o \
	calcmod.o \
	numhess.o

########################################################################
# RIXS plotting code
########################################################################
RIXSPLT = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/rixsplt/rixsmod.o \
	source/analysis/rixsplt/rixsplt.o

RIXSPLT_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	rixsmod.o \
	rixsplt.o

########################################################################
# Wavepacket autocorrelation function-to-spectrum code
########################################################################
AUTO2SPEC = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/auto2spec/auto2spec.o

AUTO2SPEC_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	auto2spec.o

########################################################################
# Chebyshev order-domain autocorrelation function-to-spectrum code
########################################################################
CHEBY2SPEC = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/cheby2spec/cheby2spec.o

CHEBY2SPEC_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	cheby2spec.o

########################################################################
# Filter diagonalisation code
########################################################################
FDIAG = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/fdiag/filtermod.o \
	source/analysis/fdiag/fdiag.o

FDIAG_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	filtermod.o \
	fdiag.o

########################################################################
# TD-NTO analysis code
########################################################################
NTOANA = $(MULTI) \
	source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/adclib/electron_density.o \
	source/iomodules/molden.o \
	source/analysis/ntoana/ntoana.o

NTOANA_OBJ = accuracy.o \
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
        fock_tools.o \
	constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	electron_density.o \
	molden.o \
	ntoana.o

########################################################################
# DPSS code
########################################################################
DPSS = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/chebyfd/dpss_ev.o \
	source/analysis/chebyfd/pythag.o \
	source/analysis/chebyfd/set_xint.o \
	source/analysis/chebyfd/sft.o \
	source/analysis/chebyfd/tinvit.o \
	source/analysis/chebyfd/tridib.o \
	source/analysis/chebyfd/xint.o \
	source/analysis/chebyfd/dpssmt.o \
	source/analysis/chebyfd/dpss.o

DPSS_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	dpss_ev.o \
	pythag.o \
	set_xint.o \
	sft.o \
	tinvit.o \
	tridib.o \
	xint.o \
	dpssmt.o \
	dpss.o 

########################################################################
# Chebyshev filter diagonalisation code
########################################################################
CHEBYFD = source/include/constants.o \
	source/include/channels.o \
	source/iomodules/iomod.o \
	source/iomodules/parsemod.o \
	source/analysis/chebyfd/dpss_ev.o \
	source/analysis/chebyfd/pythag.o \
	source/analysis/chebyfd/set_xint.o \
	source/analysis/chebyfd/sft.o \
	source/analysis/chebyfd/tinvit.o \
	source/analysis/chebyfd/tridib.o \
	source/analysis/chebyfd/xint.o \
	source/analysis/chebyfd/dpssmt.o \
	source/analysis/chebyfd/cfdmod.o \
	source/analysis/chebyfd/gaussian_coeffs.o \
	source/analysis/chebyfd/idpss_coeffs.o \
	source/analysis/chebyfd/chebyfd.o

CHEBYFD_OBJ = constants.o \
	channels.o \
	iomod.o \
	parsemod.o \
	dpss_ev.o \
	pythag.o \
	set_xint.o \
	sft.o \
	tinvit.o \
	tridib.o \
	xint.o \
	dpssmt.o \
	cfdmod.o \
	gaussian_coeffs.o \
	idpss_coeffs.o \
	chebyfd.o

#-----------------------------------------------------------------------
# Rules to create the programs
#-----------------------------------------------------------------------
adc: $(OBJECTS)
	$(F90) $(F90OPTS) $(ADC_OBJ) $(LIBS) -o bin/adc.x 
	rm -f *.o *~ *.mod 2>/dev/null

stieltjes_ap: $(STIELTJES_AP)
	$(F90) $(F90OPTS) $(STIELTJES_AP_OBJ) $(LIBS) -o bin/stieltjes_ap.x
	rm -f *.o *~ *.mod

mcspline: $(MCSPLINE)
	$(F90) $(F90OPTS) $(MCSPLINE_OBJ) -o bin/mcspline.x
	rm -f *.o *~ *.mod

numhess: $(NUMHESS)
	$(F90) $(F90OPTS) $(NUMHESS_OBJ) $(LIBS) -o bin/numhess.x
	rm -f *.o *~ *.mod

rixsplt: $(RIXSPLT)
	$(F90) $(F90OPTS) $(RIXSPLT_OBJ) $(LIBS) -o bin/rixsplt.x
	rm -f *.o *~ *.mod

auto2spec: $(AUTO2SPEC)
	$(F90) $(F90OPTS) $(AUTO2SPEC_OBJ) $(LIBS) -o bin/auto2spec.x
	rm -f *.o *~ *.mod

cheby2spec: $(CHEBY2SPEC)
	$(F90) $(F90OPTS) $(CHEBY2SPEC_OBJ) $(LIBS) -o bin/cheby2spec.x
	rm -f *.o *~ *.mod

fdiag: $(FDIAG)
	$(F90) $(F90OPTS) $(FDIAG_OBJ) $(LIBS) -o bin/fdiag.x
	rm -f *.o *~ *.mod

ntoana: $(NTOANA)
	$(F90) $(F90OPTS) $(NTOANA_OBJ) $(LIBS) -o bin/ntoana.x
	rm -f *.o *~ *.mod

dpss: $(DPSS)
	$(F90) $(F90OPTS) $(DPSS_OBJ) $(LIBS) -o bin/dpss.x
	rm -f *.o *~ *.mod

chebyfd: $(CHEBYFD)
	$(F90) $(F90OPTS) $(CHEBYFD_OBJ) $(LIBS) -o bin/chebyfd.x
	rm -f *.o *~ *.mod

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
