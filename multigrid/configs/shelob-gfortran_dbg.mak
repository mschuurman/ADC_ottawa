BUILD_ID :="Debug gfortran, built on $(shell hostname) at $(shell date)"
# corei7 with avx
CC  = gcc-4.8.2 -Wall -m64 -mavx -g -Og -march=native -mtune=native
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran-4.8.2 -I. -m64 -mavx -g -Og -march=native -mtune=native -fopenmp \
      -std=gnu -pedantic -Wall -Wno-unused-dummy-argument -Wno-unused-function \
      -Wno-maybe-uninitialized -Wcharacter-truncation -Wunderflow \
      -Wsurprising -fbacktrace -fcheck=all -fno-realloc-lhs \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
# -traceback causes crashes for some sources, which appear to be fine
F90L = $(F90) 
LAPACK = -llapack -lopenblas -lquadlapack_gfortran
FFTW3 = -lfftw3_threads -lfftw3f_threads -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/opt/dx/include
DXLIB = -L/opt/dx/lib_linux/ -lDXlite 
LIBEXTRA = # -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o

arpack_zmin.o: arpack_zmin.f
	$(F90) -c arpack_zmin.f
#	$(F90) -c -fno-check=do arpack_zmin.f
