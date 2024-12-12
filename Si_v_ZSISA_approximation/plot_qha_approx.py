#!/usr/bin/env python
"""
Quasi-harmonic approximation
============================

This example shows how to use the GSR.nc and PHDOS.nc files computed with different volumes
to compute thermodynamic properties within the quasi-harmonic approximation.
"""
import os
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.qha_aproximation import QHA_App

strains = [96,98,100,102,104,106]
strains2 = [98,100,102,104,106] #EinfVib4(D)
#strains2 = [96,98,100,102,104] #EinfVib4(S)
#strains2 = [100,102,104] #EinfVib2(D)

gsr_paths = [os.path.join( "scale_{:d}_GSR.nc".format(s)) for s in strains]
ddb_paths = [os.path.join( "scale_{:d}_GSR_DDB".format(s)) for s in strains]
dos_paths = [os.path.join( "scale_{:d}_PHDOS.nc".format(s)) for s in strains2]

#qha = QHA_App.from_files_app(gsr_paths, dos_paths)
qha = QHA_App.from_files_app_ddb(ddb_paths, dos_paths)
qha.plot_energies(tstop=800,tstart=0 ,num=11,title="Energies as a function of volume for different T")

# Vinet
qha.plot_vol_vs_t(tstop=800,tstart=0 ,num=101,title="Volume as a function of T")

qha.plot_abc_vs_t(tstop=800,tstart=0 ,num=101,title="Lattice as a function of T")
qha.plot_abc_vs_t(tstop=800,tstart=0 ,num=101,lattice="b",title="Lattice as a function of T")

qha.plot_thermal_expansion_coeff(tstop=800,tstart=0 ,num=101 ,title="Volumetric thermal expansion coefficient as a function of T")
qha.plot_thermal_expansion_coeff_abc(tstop=800,tstart=0 ,num=101 ,title="Thermal expansion coefficient as a function of T")
qha.plot_angles_vs_t(tstop=800,tstart=0, num=101,title="Angles as a function of T")

# 4th order polinomial
qha.plot_vol_vs_t_4th(tstop=800,tstart=0 ,num=101,title="Volume as a function of T")

qha.plot_abc_vs_t_4th(tstop=800,tstart=0 ,num=101,lattice="a",title="Lattice as a function of T")
qha.plot_abc_vs_t_4th(tstop=800,tstart=0 ,title="Lattice as a function of T")
qha.plot_thermal_expansion_coeff_4th( tref=293 ,title="Volumetric thermal expansion coefficient as a function of T")
qha.plot_thermal_expansion_coeff_abc_4th(tstop=800,tstart=0 ,num=101 ,tref=293,title="Thermal expansion coefficient as a function of T")
qha.plot_angles_vs_t_4th(tstop=800,tstart=0 ,num=101,angle=3,title="Angles as a function of T")
##
