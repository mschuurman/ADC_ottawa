CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = ifort -nofor_main -complex_limited_range -debug extended -traceback -warn -heap-arrays 16 \
            -ftrapuv -check all -check noarg_temp
F90L = $(F90)
LAPACK = -mkl=parallel
FFTW3 = -L/usr/lib64/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/opt/dx/include
DXLIB = -L/opt/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
