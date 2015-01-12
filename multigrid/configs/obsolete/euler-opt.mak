CC  = gcc -Wall -msse2 -mfpmath=sse -O3
F90 = ifort -warn -nofor_main -xN -axP -ipo -O3 -no-prec-div -static -openmp -assume cc_omp -complex_limited_range \
            -debug extended -traceback 
F90L = $(F90)
LAPACK = -L /opt/intel/mkl70/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lpthread
FFTW3 = -L/usr/local/lib/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/dx/include
DXLIB = -L/usr/local/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
