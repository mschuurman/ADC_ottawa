CC  = pathcc -m64
F90 = pathf90 -ansi -fullwarn -trapuv -C -m64 -O -woff1438 -woff20 -woff399
F90L = $(F90)
LAPACK = -L /opt/intel/mkl/8.1.1/lib/em64t/ -lmkl_lapack -lmkl_em64t -lguide -lpthread
FFTW3 = -L/usr/local/sims/lib64/pathscale -lfftw3 -lfftw3f
DXINC = -DUSE_DX -I/usr/local/sims/dx/include
DXLIB = -L/usr/local/sims/dx/lib_linux/ -lDXlite

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
