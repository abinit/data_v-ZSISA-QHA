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

# Initialize temperature and pressure from command-line arguments
Temp = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
Press = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
script_name = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
workdir = f"{script_name}_{Temp:04.0f}_{Press:02.0f}"

os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)

# Define strain configurations and paths
strains = [1000, 1005, 1010]

#mode='TEC'  # for thermal expansion coefficents
mode='ECs'   # for thermal expansion coefficents + elastic constants

root = os.path.join(abidata.dirpath, "data_v-ZSISA-QHA.git", "ZnO_ZSISA_ElasticConstants")
dos_paths = [[[[[[os.path.join(root, f"scale_{s1}_{s2}_{s3}_{s4}_1000_1000/out_PHDOS.nc") 
            for s6 in strains] for s5 in strains] for s4 in strains] 
            for s3 in strains] for s2 in strains] for s1 in strains]
guess_path = "Relax2o_GSR.nc"
gsr_BO_path = os.path.join(root, f"scale_1000_1000_1000_1000_1000_1000/out_GSR.nc")


def make_scf_input(qha,paral_kgb=0, stress=None, work=None):
    """
    Constructs the input file for the GS calculation.
    """
    structure = qha.structure_guess  # Take the first structure (or specify correctly)
    pseudos = ["../Zn.psp8", "../O.psp8"]
    
    gs_inp = abilab.AbinitInput(structure, pseudos=pseudos)

    # Define a dictionary for the base SCF parameters
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
        "tolvrs": 1.0e-14,  # SCF stopping criterion
    }
    # Define a dictionary for additional parameters
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

    # Set variables in gs_inp
    if work=="elastic":
        gs_inp.set_vars(**scf_vars)
    else :
        gs_inp.set_vars(**{**scf_vars, **extra_vars})

    return gs_inp
# Function to check if a specific job is running
def is_job_running(job_id):
    result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
    return job_id in result.stdout

# Function for elastic flow
def build_flow(options):
    """ Create a `Flow` for phonon calculations. """
    #options.workdir="flow_all"
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    flow = flowtk.Flow(workdir=options.workdir)

    # Build SCF input and register the first work
    scf_input = make_scf_input(qha, stress=stress, work="elastic" )
    elast_work = flowtk.ElasticWork.from_scf_input(scf_input, with_relaxed_ion=True, with_piezo=False)

    flow.register_work(elast_work)
    return flow

def compare_lattice(edata,qha):
    tolerance=1e-6        
    if np.allclose(edata.structure.lattice.matrix, qha.structure_guess.lattice.matrix, atol=tolerance):
        lat_match=True
    else:
        print("The lattice matrices are different.")

    for site1, site2 in zip(edata.structure.sites, qha.structure_guess.sites):
        if site1.species != site2.species:
            print("Mismatch in atomic species:", site1.species, site2.species)
        elif not np.allclose(site1.frac_coords, site2.frac_coords, atol=tolerance):
            print("Mismatch in fractional coordinates:", site1.frac_coords, site2.frac_coords)
        elif not np.allclose(site1.coords, site2.coords, atol=tolerance):
            print("Mismatch in cartesian coordinates:", site1.coords, site2.coords)
        else:
            site_match=True
    if lat_match and site_match :
        check=True 
    return check

while True:
    try:
        print(f"Calculating stress for T={Temp} and pressure={Press}...")

        # Initialize QHA_ZSISA object and compute stress
        qha = QHA_ZSISA.from_files(dos_paths, guess_path, gsr_BO_path)
        result, stress = qha.cal_stress(Temp, Press, mode)
        # Generate SCF input
        gs_inp = make_scf_input(qha, stress=stress )
        gs_inp.write("Relax2.abi")
        prefix_to_remove = "prefix"
        subprocess.run(f"sed -i '/{prefix_to_remove}/d' Relax2.abi", shell=True)

        if result:
            print("Stress calculation complete. Initiating elastic calculations...")
            workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")
            ddbdir=workdir+"/w0/outdata/out_DDB"
            if  os.path.exists(ddbdir):
                with DdbFile(ddbdir) as ddb:
                    edata = ddb.anaget_elastic()
                    check=compare_lattice(edata,qha)
                    if check:
                        with open("elastic_constant.txt", "w") as f:
                            f.write(str(edata))
            else:
                shutil.rmtree(workdir)
                        
            if not os.path.exists(ddbdir):
                # Parse command-line options and run phonon calculations
                options = flowtk.build_flow_main_parser().parse_args([])
                flow = build_flow(options)
                flow.make_scheduler().start()
                with DdbFile(ddbdir) as ddb:
                    edata = ddb.anaget_elastic()
                    with open("elastic_constant.txt", "w") as f:
                        f.write(str(edata))

            result, stress = qha.cal_stress(Temp, Press, mode )
            break  # Exit the loop

        # Submit the job using sbatch
        sbatch_result = subprocess.run(["sbatch", "../job.sh"], capture_output=True, text=True)
        job_id = sbatch_result.stdout.strip().split()[-1]  # Get the job ID from sbatch output
        print("Submitted batch job", job_id)

        # Wait for the job to finish
        while is_job_running(job_id):
            print("Job", job_id, "is still running...")
            time.sleep(30)  # Check every 30 seconds

        print("Job", job_id, "has completed.")

    except subprocess.CalledProcessError as e:
        print("Error during execution:", e)

    print("Restarting the process...")

