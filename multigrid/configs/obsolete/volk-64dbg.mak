CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = ifort -nofor_main -C -debug extended -traceback -check noarg_temp_created \
            -fpstkchk -ftrapuv -warn -error_limit 5 -xP -O2 -complex_limited_range
#           -openmp -assume cc_omp 
# -traceback causes crashes for some sources, which appear to be fine
F90L = $(F90)
LAPACK = -L /opt/intel/mkl/9.0/lib/em64t/ -lmkl_lapack -lmkl_em64t -lguide -lpthread
FFTW3 = -L/usr/local/lib64/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/dx/include
DXLIB = -L/usr/local/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
