#F90	=/opt/intel/bin/ifort
#CC	=/opt/intel/bin/icc
F90	= gfortran
F77	= gfortran
CC	= gcc

#F90OPTS =  -g -CB 
#CCOPTS  =  -g -O0

F90OPTS = -g -ffixed-line-length-none -ffree-line-length-none
CCOPTS  = -g -O0

# export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
#  export LD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64:$LD_LIBRARY_PATH
#LIBS = -L/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64  /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_blas95_lp64.a  /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm   /home/mr2911/phys/phisok/phis/lib/libphis.a -L/home/mr2911/Molcas/molcas76/lib -lmolcas -L/home/mr2911/Molcas/molcas76/g/lib/LINUX64 -lma -L/home/mr2911/Francesco.LIB -lblnz -ldinvop -lmem -lutil #-L/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64  /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_blas95_lp64.a  /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm   #/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_intel_lp64.a  /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_core.a /opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64/libmkl_sequential.a  -L/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64 -lpthread  -lm
# /cvos/shared/TC/phis/molcas/74.ifc/lib/libphis.a -L/cvos/shared/TC/molcas/serial/molcas74.ifc/lib -lmolcas -L/cvos/shared/TC/molcas/serial/molcas74.ifc/g/lib/LINUX64 -lma \
# -L/home/soeren/libs/franc_libs -lblnz -ldinvop -lmem -lutil -lguide -L/cvos/shared/apps/intel/mkl/10.1.3.027/lib/em64t -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread
MDIR=multigrid_interface

LIBS= -L/usr/lib64 -lblas -llapack

#MULTI=$(MDIR)/accuracy.o $(MDIR)/timer.o $(MDIR)/lapack.o $(MDIR)/dgefa.o $(MDIR)/dgedi.o $(MDIR)/math.o $(MDIR)/matrix_tools.o $(MDIR)/gamess_internal.o $(MDIR)/import_gamess.o $(MDIR)/os_integral_operators.o $(MDIR)/integral_tools.o $(MDIR)/integrals_mo2e.o

MULTI=accuracy.o printing.o timer.o lapack.o dgefa.o dgedi.o math.o matrix_tools.o os_integral_operators.o gamess_internal.o import_gamess.o integral_tools.o integrals_mo2e.o
ADC = constants.o parameters.o misc.o external.o external_diag.o filetools.o adc_ph.o  dipole_ph.o D_matrix.o read_param.o sym_allowed_exc.o select_fano.o get_matrix.o  get_matrix_dipole_complete.o get_moment.o fspacetrial.o fspace2.o davmod.o photoionisation.o   Propagate.o   master_adc1_prop.o    master_adc2_prop.o   master_adc2ext_prop.o  main_draft1.o

OBJECTS=$(MULTI) $(ADC)

ww: $(OBJECTS)
	$(F90) $(F90OPTS) $(OBJECTS) $(LIBS) -o  adc.x 

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean:
	rm -f *.o *~ *.mod 
