#!/usr/bin/env python
"""
Quasi-harmonic approximation
============================

This example shows how to use the GSR.nc and PHDOS.nc files computed with different volumes
to compute thermodynamic properties within the quasi-harmonic approximation.
"""
import os
import numpy as np
import abipy.data as abidata

from abipy.dfpt.phonons import PhononBands
from abipy.dfpt.qha_2D import QHA_2D

bo_strains_a = [995, 1000, 1005, 1010, 1015]
bo_strains_c = [995, 1000, 1005, 1010, 1015]
phdos_strains_a = [1000, 1005, 1010]
phdos_strains_c = [1000, 1005, 1010]

#gsr_paths = [[f"scale_{s1}_{s3}/out_GSR.nc" for s3 in bo_strains_c] for s1 in bo_strains_a]
gsr_paths = [[f"scale_{s1}_{s3}/out_GSR_DDB" for s3 in bo_strains_c] for s1 in bo_strains_a]
phdos_paths = [[f"scale_{s1}_{s3}/out_PHDOS.nc" for s3 in phdos_strains_c] for s1 in phdos_strains_a]

bo_strains_ac = [bo_strains_a, bo_strains_c]
bo_strains_ac = (np.array(bo_strains_ac) - 1000) / 100
phdos_strains_ac = [phdos_strains_a, phdos_strains_c]
phdos_strains_ac = (np.array(phdos_strains_ac) - 1000) / 100

#qha = QHA_2D.from_files(gsr_paths, phdos_paths, bo_strains_ac, phdos_strains_ac, gsr_file="GSR.nc")
qha = QHA_2D.from_files(gsr_paths, phdos_paths, bo_strains_ac, phdos_strains_ac, gsr_file="DDB")

qha.plot_energies()
qha.plot_free_energies(tstop=500, tstart=0, num=6)
qha.plot_thermal_expansion(tstop=1000, tstart=0, num=101)
qha.plot_lattice(tstop=1000, tstart=00, num=101)
