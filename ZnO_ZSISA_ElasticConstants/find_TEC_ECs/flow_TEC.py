# === Import required modules ===
import subprocess
import time
import sys
import os
import numpy as np
import shutil
import abipy.data as abidata

from abipy.dfpt.qha_general_stress import QHA_ZSISA
from abipy import abilab
from abipy import flowtk
from abipy.dfpt.ddb import DdbFile

"""
Script for Iterative Determination of Lattice Parameters at Finite Temperature and External Pressure.

This script implements the iterative process described in the flowchart for obtaining the optimal
lattice configuration at a given temperature (T) and external pressure (Pext). The procedure is
based on first-principles calculations and the ZSISA-E\infty vib2 framework, integrating both thermal
and Born–Oppenheimer (BO) contributions to stress.

Workflow Summary:
1. An initial guess for the lattice configuration [R] is provided.
2. Phonon density of states (DOS) and static DFT results are used to compute the thermal stress
   at the target temperature and external pressure.
3. A target stress is constructed from the thermal stress and Pext.
4. A DFT structural relaxation (with cell degrees of freedom) is performed using this target stress.
5. The new relaxed structure is used to recompute the stress tensor.
6. Steps 2–5 are repeated iteratively until the BO stress matches the target stress within tolerance.
7. Once convergence is achieved, a final elastic constants calculation is triggered using the
   converged relaxed structure.

The script automatically handles job submission, monitoring, and result extraction for
relaxation and elasticity calculations using ABINIT, via AbiPy's flow toolkit.
"""

# === Read temperature and pressure from command-line arguments ===
Temp = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
Press = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

# === Construct working directory name based on script name and inputs ===
script_name = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
workdir = f"{script_name}_{Temp:04.0f}_{Press:02.0f}"

# === Create and enter the working directory ===
os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)

# === Define strain scaling factors for QHA grid ===
strains = [1000, 1005, 1010]

# === QHA mode: 'TEC' for thermal expansion only, 'ECs' to include elastic constants ===
mode = 'ECs'

# === Define paths to phonon DOS files for all strain grid points ===
root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_ElasticConstants")
dos_paths = [[[[[[os.path.join(root, f"scale_{s1}_{s2}_{s3}_{s4}_1000_1000/out_PHDOS.nc") 
            for s6 in strains] for s5 in strains] for s4 in strains] 
            for s3 in strains] for s2 in strains] for s1 in strains]

# === Paths to guessed and BO-relaxed structures ===
guess_path = "Relax2o_GSR.nc"
gsr_BO_path = os.path.join(root, f"scale_1000_1000_1000_1000_1000_1000/out_GSR.nc")

# === Function to generate SCF input file with or without target stress ===
def make_scf_input(qha, paral_kgb=0, stress=None, work=None):
    structure = qha.structure_guess
    pseudos = ["../Zn.psp8", "../O.psp8"]

    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    # SCF settings
    scf_vars = {
        "paral_kgb": paral_kgb,
        "nline": 10,
        "nbdbuf": 0,
        "nshiftk": 1,
        "nspden": 1,
        "ngkpt": [6, 6, 4],
        "shiftk": [0.0, 0.0, 0.0],
        "charge": 0.0,
        "nstep": 100,
        "ecut": 42.0,
        "ecutsm": 1.0,
        "chksymbreak": 1,
        "occopt": 1,
        "nband": 26,
        "kptopt": 1,
        "tolvrs": 1.0e-14,
    }

    # Relaxation settings (only used if work is not "elastic")
    extra_vars = {
        "strtarget": stress,
        "ionmov": 2,
        "ntime": 60,
        "optcell": 2,
        "dilatmx": 1.04,
        "tolmxf": 1.0e-5,
        "strfact": 1000.,
        "prtden": 0,
        "prtwf": 0,
        "prteig": 0
    }

    # Merge variables based on whether it's a relaxation or elasticity run
    if work == "elastic":
        gs_inp.set_vars(**scf_vars)
    else:
        gs_inp.set_vars(**{**scf_vars, **extra_vars})

    return gs_inp

# === Check whether a SLURM job is still running ===
def is_job_running(job_id):
    result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
    return job_id in result.stdout

# === Build elastic constants workflow ===
def build_flow(options):
    flow = flowtk.Flow(workdir="elastic")
    scf_input = make_scf_input(qha, stress=stress, work="elastic")
    elast_work = flowtk.ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=False)
    flow.register_work(elast_work)
    return flow

# === Compare relaxed structure from DDB with QHA guess ===
def compare_lattice(edata, qha):
    tolerance = 1e-6
    lat_match = np.allclose(edata.structure.lattice.matrix, qha.structure_guess.lattice.matrix, atol=tolerance)
    site_match = True

    for site1, site2 in zip(edata.structure.sites, qha.structure_guess.sites):
        if site1.species != site2.species:
            print("Mismatch in atomic species:", site1.species, site2.species)
            site_match = False
        elif not np.allclose(site1.frac_coords, site2.frac_coords, atol=tolerance):
            print("Mismatch in fractional coordinates:", site1.frac_coords, site2.frac_coords)
            site_match = False
        elif not np.allclose(site1.coords, site2.coords, atol=tolerance):
            print("Mismatch in cartesian coordinates:", site1.coords, site2.coords)
            site_match = False

    return lat_match and site_match

# === Main loop to calculate stress and trigger elastic constants ===
while True:
    try:
        print(f"Calculating stress for T={Temp} and pressure={Press}...")

        # Initialize QHA object and compute target stress
        qha = QHA_ZSISA.from_files(dos_paths, guess_path, gsr_BO_path)
        chk_converge, stress = qha.cal_stress(Temp, Press, mode)

        # Write input for relaxation
        gs_inp = make_scf_input(qha, stress=stress)
        gs_inp.write("Relax2.abi")
        subprocess.run(f"sed -i '/prefix/d' Relax2.abi", shell=True)

        if chk_converge:
            print("Stress calculation complete. Initiating elastic calculations...")

            # Set working dir and expected DDB path
            workdir = "elastic"
            ddbdir = workdir + "/w0/outdata/out_DDB"

            # Check if DDB already exists and is consistent
            if os.path.exists(workdir):
                if os.path.exists(ddbdir):
                    with DdbFile(ddbdir) as ddb:
                        edata = ddb.anaget_elastic()
                        if compare_lattice(edata, qha):
                            with open("elastic_constant.txt", "w") as f:
                                f.write(str(edata))
                else:
                    shutil.rmtree(workdir)  # Clean incomplete output

            # Run elastic constants workflow if needed
            if not os.path.exists(ddbdir):
                options = flowtk.build_flow_main_parser().parse_args([])
                flow = build_flow(options)
                flow.make_scheduler().start()

                # Wait for DDB and extract elastic constants
                with DdbFile(ddbdir) as ddb:
                    edata = ddb.anaget_elastic()
                    with open("elastic_constant.txt", "w") as f:
                        f.write(str(edata))

            # Recalculate stress to confirm convergence
            chk_converge, stress = qha.cal_stress(Temp, Press, mode)
            break  # Done!

        # === Submit relaxation job if not converged ===
        sbatch_result = subprocess.run(["sbatch", "../job.sh"], capture_output=True, text=True)
        job_id = sbatch_result.stdout.strip().split()[-1]
        print("Submitted batch job", job_id)

        # Monitor job until completion
        while is_job_running(job_id):
            print("Job", job_id, "is still running...")
            time.sleep(30)

        print("Job", job_id, "has completed.")

    except subprocess.CalledProcessError as e:
        print("Error during execution:", e)

    print("Restarting the process...")

