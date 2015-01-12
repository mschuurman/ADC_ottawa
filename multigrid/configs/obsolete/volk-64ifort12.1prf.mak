CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = /opt/intel/composer_xe_2011_sp1.8.273/bin/intel64/ifort \
            -warn -nofor_main -assume buffered_io -xSSE3 -axSSE4.2 -complex_limited_range \
            -O3 -ipo -no-prec-div -static -fp-model fast -heap-arrays 32 -pad \
            -mkl=sequential -opt-prefetch \
            -debug extended -traceback -p
F90L = $(F90)
LAPACK = 
FFTW3 = -L/usr/local/lib64/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/dx/include
DXLIB = -L/usr/local/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
