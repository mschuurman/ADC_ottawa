# corei7 or barcelona
CC  = gcc -Wall -m64 -msse2 -g -O0 -march=native -mtune=native
F90 = gfortran -m64 -msse2 -g -O0 -march=native -mtune=native -fopenmp \
      -std=gnu -pedantic -Wall -Wno-unused-dummy-argument -Wcharacter-truncation -Wunderflow \
      -fbacktrace -fcheck=all -fno-realloc-lhs
# -traceback causes crashes for some sources, which appear to be fine
F90L = $(F90) 
LAPACK = /usr/lib64/gcc/x86_64-suse-linux/4.7/libgomp.a \
         -L /opt/intel/mkl/9.0/lib/em64t/ -lmkl_lapack -lmkl_em64t -lguide -lpthread -lrt -lc \
	 -L /opt/intel/fce/9.1.041/lib/ -lsvml -lirc
FFTW3 = -L/usr/local/lib64/ -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/dx/include
DXLIB = -L/usr/local/dx/lib_linux/ -lDXlite 

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c -fno-openmp fftw_fftw3.f90 -o fftw.o

arpack_zmin.o: arpack_zmin.f
	$(F90) -c arpack_zmin.f
#	$(F90) -c -fno-check=do arpack_zmin.f
