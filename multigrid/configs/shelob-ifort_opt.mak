BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
CC  = gcc -Wall -m64 -march=k8 -mavx -mfpmath=sse -O3
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = /opt/intel/composer_xe_2013_sp1.3.174/bin/intel64/ifort \
            -static -warn -nofor_main -assume buffered_io \
            -O3 -ip -no-prec-div -xAVX -complex_limited_range -fp-model fast -openmp -heap-arrays 32 -pad \
            -mkl=sequential -openmp-link=static \
            -opt-matmul -no-opt-prefetch -finline-limit=500 \
            -debug extended -traceback \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90)
LAPACK = -L/home/ps/lib64 -lquadlapack_intel
FFTW3 = -L/home/ps/lib/intel/lib64 -Wl,--whole-archive -lfftw3 -lfftw3f -Wl,--no-whole-archive -lfftw3_threads -lfftw3f_threads
DXINC = -DUSE_DX -I/opt/dx/include
DXLIB = -L/opt/dx/lib_linux/ -lDXlite
LIBEXTRA = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o

arpack_zmin.o: arpack_zmin.f
	$(F90) -c arpack_zmin.f

