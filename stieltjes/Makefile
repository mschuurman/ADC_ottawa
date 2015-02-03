########################################################################
#
#                  MAKEFILE FOR STIELTJES IMAGING CODE                  
#
########################################################################

#-----------------------------------------------------------------------
# Compiler flags
#
# N.B. we now have to use the C preprocessor due to interfacing with
# PETSc/SLEPc
#-----------------------------------------------------------------------
F90	= gfortran
F77	= gfortran
CC	= gcc

F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none
CCOPTS  = -g -O0

#-----------------------------------------------------------------------
# Define object files
#-----------------------------------------------------------------------
OBJECTS = qmath.o \
	pythag_quad.o \
	tql2_quad.o \
        globalmod.o \
	stieltjes_s3_modified.o \
	main_stieltjes.o 

#-----------------------------------------------------------------------
# Rules to create the program
#-----------------------------------------------------------------------
stieltjes: $(OBJECTS)
	$(F90) $(F90OPTS) $(OBJECTS) -o stieltjes.x

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod 
