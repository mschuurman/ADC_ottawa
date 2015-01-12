# corei7 or barcelona
CC  = gcc -Wall -m64 -msse2 -g -O0 -march=native -mtune=native
F90 = gfortran -m64 -msse2 -g -O0 -march=native -mtune=native \
      -std=gnu -pedantic -Wall -Wno-unused-dummy-argument -Wcharacter-truncation -Wunderflow \
      -fbacktrace -fcheck=all -fno-realloc-lhs
# -traceback causes crashes for some sources, which appear to be fine
F90L = $(F90) 
LAPACK = -llapack -lblas
FFTW3 = -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/opt/dx/include
DXLIB = -L/opt/dx/lib_linux/ -lDXlite 

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
