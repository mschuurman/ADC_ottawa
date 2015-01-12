goal: makefile.dep
	# make crude_adiabatic_ff.x
	# make scf_ca.x
	# make test_2e.x
	# make partial_waves_sp.x
	# make mo_adk.x
	# make partial_waves_c3.x
	# make amplitudes.x
	# make amplitudes_eva.x
	# make eikonal_integrals.x
	# make eikonal_integrals_2index.x
	# make h2p_recollide_3d.x
	# make partition-in-small-box.x
	# make ee_integral.x
	# make ecp_test.x
	# make correlation_potential_v2.x
	# make correlation_potential_v3.x
	# make correlation_potential.x
	# make dyson_nchan.x
	# make dyson_1chan.x
	# make dyson_harmonics.x
	# make spherical_interpolate.x
	# make sample_orbital.x
	# make chiral.x
	# make orbital_fft.x
	# make orbital_fft_v2.x
	# make symgen.x
	 make static_tunnel.x
	# make basis_recontract.x

MAKEFLAGS = -r

.SUFFIXES: .f90 .o .x .c .dep

#
# Default settings
#
EXTNAMES = -Ddriver=driver_ -DdrvGridHead=drvgridhead_ -DdrvGridData=drvgriddata_ \
	   -DdrvGridShow=drvgridshow_

LIBARPACK = # -larpack
OBJARPACK = arpack_zmin.o

arpack_zmin.o: arpack_zmin.f
	$(F90) -c arpack_zmin.f

#
#  This is the default; a config file may override it.
#
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
#ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements

#
# System-specific overrides
#
# include configs/shelob-ifort_opt.mak
# include configs/shelob-gfortran_opt.mak
# include configs/shelob-gfortran_dbg.mak
  include configs/shelob-gfortran_noquad_opt.mak
# include configs/shelob-gfortran_noquad_dbg.mak

#
# Finish the set-up
#
LIBS = $(LIBARPACK) $(LAPACK) $(FFTW3) $(DXLIB) $(LAPACK) $(LIBEXTRA)

#
# Compiling and archiving rules
#
.f90.o:
	$(ACT) $< >preprocess/$<
	$(F90) -c preprocess/$<

dgefa.o:        dgefa.f
	$(F90) -c dgefa.f

dgedi.o:        dgedi.f
	$(F90) -c dgedi.f

erf.o:          erf.f
	$(F90) -c erf.f

derf.o:         derf.f
	$(F90) -c derf.f

struve.o:       struve.f
	$(F90) -c struve.f

.c.o:
	$(CC) -c $<

clean:
	-/bin/rm -f *.{o,mod,x,il,a} *__genmod.f90 checkpoint_{field,main}.* makefile.dep

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Interactive or batch online interface
#
#driver.o:	onlinedriver.c
#	$(CC) $(CFLAGS) $(DXINC) $(EXTNAMES) -c onlinedriver.c -o driver.o

driver.o:	exportdriver.c
	$(CC) $(CFLAGS) $(DXINC) $(EXTNAMES) -c exportdriver.c -o driver.o

#
# Explicit dependencies
#

LIBMULTIMIN += $(OBJARPACK)
LIBMULTIMIN += accuracy.o
LIBMULTIMIN += arpack.o
LIBMULTIMIN += atoms.o
LIBMULTIMIN += basis_cap.o
LIBMULTIMIN += bessel.o
LIBMULTIMIN += biorthogonal_tools.o
LIBMULTIMIN += block_diag.o
LIBMULTIMIN += ci_tools.o
LIBMULTIMIN += complex.o
LIBMULTIMIN += cubeint.o
LIBMULTIMIN += diis.o
LIBMULTIMIN += dft.o
LIBMULTIMIN += dgedi.o
LIBMULTIMIN += dgefa.o
LIBMULTIMIN += ecp_convert_gamess.o
LIBMULTIMIN += fftw.o
LIBMULTIMIN += fields.o
LIBMULTIMIN += fock_tools.o
LIBMULTIMIN += gamess_internal.o
LIBMULTIMIN += import_gamess.o
LIBMULTIMIN += integral_tools.o
LIBMULTIMIN += integrals_mo2e.o
LIBMULTIMIN += interpolate.o
LIBMULTIMIN += lapack.o
LIBMULTIMIN += lebedev.o
LIBMULTIMIN += math.o
LIBMULTIMIN += matrix_tools.o
LIBMULTIMIN += molecular_grid.o
LIBMULTIMIN += mp2_tools.o
LIBMULTIMIN += opendx.o
LIBMULTIMIN += os_integral_operators.o
LIBMULTIMIN += parallel.o
LIBMULTIMIN += point_group.o
LIBMULTIMIN += pulse.o
LIBMULTIMIN += printing.o
LIBMULTIMIN += rf_cap.o
LIBMULTIMIN += scf_tools.o
LIBMULTIMIN += sfa.o
LIBMULTIMIN += sort_tools.o
LIBMULTIMIN += struve.o
LIBMULTIMIN += timer.o

LIBMULTIGRID = $(LIBMULTIMIN)
LIBMULTIGRID += caps.o
LIBMULTIGRID += dyson_tools.o
LIBMULTIGRID += eikonal_tools.o
LIBMULTIGRID += eikonal_tools_eva.o
LIBMULTIGRID += ecp_gamess.o
LIBMULTIGRID += flux.o
LIBMULTIGRID += fock.o
LIBMULTIGRID += kinetic.o
LIBMULTIGRID += lanczos_complex.o
LIBMULTIGRID += lanczos.o
LIBMULTIGRID += liu.o
LIBMULTIGRID += multigrid.o
LIBMULTIGRID += phase_interpolation.o
LIBMULTIGRID += phase_interpolation_eva.o
LIBMULTIGRID += qmech.o
LIBMULTIGRID += symmetry.o
LIBMULTIGRID += vandalyse.o

LIBPARTITION += accuracy.o
LIBPARTITION += convolution.o
LIBPARTITION += dgedi.o
LIBPARTITION += dgefa.o
LIBPARTITION += h2_thermo.o
LIBPARTITION += h2_potential.o
LIBPARTITION += h2_rotation.o
LIBPARTITION += fftw.o
LIBPARTITION += lapack.o
LIBPARTITION += math.o
LIBPARTITION += printing.o
LIBPARTITION += sparse.o
LIBPARTITION += timer.o
LIBPARTITION += vector.o

#
# Building the binaries
#
h2p_recollide_3d.x: h2p_recollide_3d.o driver.o $(LIBMULTIGRID)
	$(F90) -o h2p_recollide_3d.x driver.o h2p_recollide_3d.o $(LIBMULTIGRID) $(LIBS)

partial_waves_sp.x: partial_waves_sp.o driver.o $(LIBMULTIGRID)
	$(F90) -o partial_waves_sp.x driver.o partial_waves_sp.o $(LIBMULTIGRID) $(LIBS)

chiral.x: chiral.o $(LIBMULTIMIN)
	$(F90) -o chiral.x chiral.o $(LIBMULTIMIN) $(LIBS)

symgen.x: symgen.o $(LIBMULTIMIN)
	$(F90) -o symgen.x symgen.o $(LIBMULTIMIN) $(LIBS)

mo_adk.x: mo_adk.o driver.o $(LIBMULTIGRID)
	$(F90) -o mo_adk.x driver.o mo_adk.o $(LIBMULTIGRID) $(LIBS)

partial_waves_c3.x: partial_waves_c3.o driver.o $(LIBMULTIGRID)
	$(F90) -o partial_waves_c3.x driver.o partial_waves_c3.o $(LIBMULTIGRID) $(LIBS)

eikonal_integrals.x: eikonal_integrals.o driver.o $(LIBMULTIGRID)
	$(F90) -o eikonal_integrals.x driver.o eikonal_integrals.o $(LIBMULTIGRID) $(LIBS)

eikonal_integrals_2index.x: eikonal_integrals_2index.o driver.o $(LIBMULTIGRID)
	$(F90) -o eikonal_integrals_2index.x driver.o eikonal_integrals_2index.o $(LIBMULTIGRID) $(LIBS)

orbital_fft.x: orbital_fft.o driver.o $(LIBMULTIGRID)
	$(F90) -o orbital_fft.x driver.o orbital_fft.o $(LIBMULTIGRID) $(LIBS)

orbital_fft_v2.x: orbital_fft_v2.o driver.o $(LIBMULTIGRID)
	$(F90) -o orbital_fft_v2.x driver.o orbital_fft_v2.o $(LIBMULTIGRID) $(LIBS)

amplitudes.x: amplitudes.o driver.o $(LIBMULTIGRID)
	$(F90) -o amplitudes.x driver.o amplitudes.o $(LIBMULTIGRID) $(LIBS)

dyson_1chan.x: dyson_1chan.o driver.o $(LIBMULTIGRID)
	$(F90) -o dyson_1chan.x driver.o dyson_1chan.o $(LIBMULTIGRID) $(LIBS)

dyson_nchan.x: dyson_nchan.o driver.o $(LIBMULTIGRID)
	$(F90) -o dyson_nchan.x driver.o dyson_nchan.o $(LIBMULTIGRID) $(LIBS)

ecp_test.x: ecp_test.o driver.o $(LIBMULTIGRID)
	$(F90) -o ecp_test.x driver.o ecp_test.o $(LIBMULTIGRID) $(LIBS)

test_2e.x: test_2e.o driver.o $(LIBMULTIGRID)
	$(F90) -o test_2e.x driver.o test_2e.o $(LIBMULTIGRID) $(LIBS)

scf_ca.x: scf_ca.o driver.o $(LIBMULTIMIN)
	$(F90) -o scf_ca.x scf_ca.o $(LIBMULTIMIN) $(LIBS)

static_tunnel.x: static_tunnel.o $(LIBMULTIMIN)
	$(F90) -o static_tunnel.x static_tunnel.o $(LIBMULTIMIN) $(LIBS)

basis_recontract.x: basis_recontract.o $(LIBMULTIMIN)
	$(F90) -o basis_recontract.x basis_recontract.o $(LIBMULTIMIN) $(LIBS)

ee_integral.x: ee_integral.o driver.o $(LIBMULTIGRID)
	$(F90) -o ee_integral.x driver.o ee_integral.o $(LIBMULTIGRID) $(LIBS)

dyson_harmonics.x: dyson_harmonics.o driver.o $(LIBMULTIGRID)
	$(F90) -o dyson_harmonics.x driver.o dyson_harmonics.o $(LIBMULTIGRID) $(LIBS)

spherical_interpolate.x: spherical_interpolate.o driver.o $(LIBMULTIGRID)
	$(F90) -o spherical_interpolate.x driver.o spherical_interpolate.o $(LIBMULTIGRID) $(LIBS)

amplitudes_eva.x: amplitudes_eva.o driver.o $(LIBMULTIGRID)
	$(F90) -o amplitudes_eva.x driver.o amplitudes_eva.o $(LIBMULTIGRID) $(LIBS)

correlation_potential.x: correlation_potential.o driver.o $(LIBMULTIGRID)
	$(F90) -o correlation_potential.x driver.o correlation_potential.o $(LIBMULTIGRID) $(LIBS)

correlation_potential_v2.x: correlation_potential_v2.o driver.o $(LIBMULTIGRID)
	$(F90) -o correlation_potential_v2.x driver.o correlation_potential_v2.o $(LIBMULTIGRID) $(LIBS)

correlation_potential_v3.x: correlation_potential_v3.o driver.o $(LIBMULTIGRID)
	$(F90) -o correlation_potential_v3.x driver.o correlation_potential_v3.o $(LIBMULTIGRID) $(LIBS)

crude_adiabatic_ff.x: crude_adiabatic_ff.o $(LIBMULTIMIN)
	$(F90) -o crude_adiabatic_ff.x crude_adiabatic_ff.o $(LIBMULTIMIN) $(LIBS)

sample_orbital.x: sample_orbital.o driver.o $(LIBMULTIGRID)
	$(F90) -o sample_orbital.x driver.o sample_orbital.o $(LIBMULTIGRID) $(LIBS)

partition-in-small-box.x: partition.o $(LIBPARTITION) 
	$(F90) -o partition-in-small-box.x partition.o $(LIBPARTITION) $(LIBS)

#
# Automatically-generated dependencies
#
include makefile.dep

