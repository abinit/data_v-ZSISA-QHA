import subprocess
import time
import sys
import os

# Import QHA_ZSISA class from abipy
from abipy.dfpt.qha_general_stress import QHA_ZSISA

# Function to check if a specific job is running
def is_job_running(job_id):
    result = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
    return job_id in result.stdout

# Initialize temperature and pressure from command-line arguments
Temp = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
Press = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

# Write initial iteration state to file
with open("iter_file.txt", "w") as f:
    f.write("0")

# Define strain configurations and paths
strains_a = [1000, 1005, 1010]
strains_c = [1000, 1005, 1010]

gsr_paths = [[f"../../scale_{s1}_{s3}/out_GSR_DDB" for s3 in strains_c] for s1 in strains_a]
dos_paths = [[f"../../scale_{s1}_{s3}/out_PHDOS.nc"     for s3 in strains_c] for s1 in strains_a]
guess_path = ["Relax2o_GSR.nc"]


# Main processing loop
while True:
    try:
        print(f"Calculating stress for T={Temp} and pressure={Press}...")

        # Initialize QHA_ZSISA object
        qha = QHA_ZSISA.from_files(gsr_paths, dos_paths, guess_path)
        result = qha.cal_stress(Temp, Press)
        print("Stress calculation result:", result)

        # Prepare Relax2.abi file
        bash_command = "cat ../temp_relax.abi lat_info.txt > Relax2.abi"
        subprocess.run(bash_command, shell=True)


        # Submit the job using sbatch
        sbatch_result = subprocess.run(["sbatch", "../job.sh"], capture_output=True, text=True)
        job_id = sbatch_result.stdout.strip().split()[-1]  # Get the job ID from sbatch output
        print("Submitted batch job", job_id)

        # Wait for the job to finish
        while is_job_running(job_id):
            print("Job", job_id, "is still running...")
            time.sleep(30)  # Check every 30 seconds

        print("Job", job_id, "has completed.")
        if result:
            # Save results with temperature-specific suffix
            bash_command0 = f"cp Relax2.abi Relax2.abi{Temp}"
            bash_command1 = f"cp Relax2o_DDB Relax2o_DDB{Temp}"
            subprocess.run(bash_command0, shell=True)
            subprocess.run(bash_command1, shell=True)

            with open("iter_file.txt", "w") as f:
                f.write("1")  # Mark as completed
            print("Calculation completed successfully. Exiting loop.")
            break  # Exit the loop

    except subprocess.CalledProcessError as e:
        print("Error during execution:", e)

    print("Restarting the process...")

