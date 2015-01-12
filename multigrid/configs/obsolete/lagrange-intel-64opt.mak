CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = ifort -warn -nofor_main -fast -openmp -assume cc_omp -complex_limited_range \
            -debug extended -traceback 
F90L = $(F90)
LAPACK = -L /opt/intel/mkl/8.1.1/lib/em64t/ -lmkl_lapack -lmkl_em64t -lguide -lpthread
FFTW3 = -L/usr/local/sims/lib64/pathscale -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/sims/dx/include
DXLIB = -L/usr/local/sims/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
