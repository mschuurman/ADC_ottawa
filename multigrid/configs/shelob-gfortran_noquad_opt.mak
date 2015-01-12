BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
# corei7 with avx
CC  = gcc -Wall -m64 -mavx -O3 -march=native -mtune=native
# -floop-interchange and -floop-strip-mine seem to cause an
# unacceptable increase in numerical errors; since they also
# give a marginal to non-existent speedup, we won't use them
#     -floop-interchange -floop-strip-mine
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
F90 = gfortran -I. -m64 -mavx -O3 -march=native -mtune=native -fopenmp \
      -floop-block \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -static -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
# -traceback causes crashes for some sources, which appear to be fine
F90L = $(F90) 
LAPACK = -llapack -lopenblas 
FFTW3 = -lfftw3_threads -lfftw3f_threads -lfftw3 -lfftw3f
#DXINC = -DUSE_DX -I/opt/dx/include
#DXLIB = -L/opt/dx/lib_linux/ -lDXlite 
LIBEXTRA = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o

arpack_zmin.o: arpack_zmin.f
	$(F90) -c arpack_zmin.f
#	$(F90) -c -fno-check=do arpack_zmin.f
