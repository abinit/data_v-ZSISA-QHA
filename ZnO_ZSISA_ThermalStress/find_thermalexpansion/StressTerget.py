#!/usr/bin/env python
import os
import sys
import subprocess
from abipy.dfpt.qha_general_stress import QHA_ZSISA


with open("iter_file.txt", "w") as f:
    f.write (f"{0}")
strains_a = [1000, 1005, 1010]
strains_c = [1000, 1005, 1010]

gsr_paths = [[f"../../scale_{s1}_{s3}/out_GSR_DDB"    for s3 in strains_c]  for s1 in strains_a]
dos_paths = [[f"../../scale_{s1}_{s3}/out_PHDOS.nc"   for s3 in strains_c]  for s1 in strains_a]

guess_path = ["Relax2o_GSR.nc" ]

qha = QHA_ZSISA.from_files(gsr_paths, dos_paths , guess_path)
pressure=0.0 
Temp=0.0 
if (len(sys.argv)>1):
    Temp=float(sys.argv[1]) 
if (len(sys.argv)>2):
    pressure=float(sys.argv[2]) 

file0="lat_info.txt" + str(Temp)
#print("T=",Temp)
result = qha.cal_stress(Temp , pressure)
print (result)

bash_command = "cat ../temp_relax.abi lat_info.txt > Relax2.abi"
subprocess.run(bash_command, shell=True)

if result==True :
    bash_command0 =  "cp Relax2.abi Relax2.abi" + str(Temp)
    bash_command2 =  "cp lat_info.txt lat_info.txt" + str(Temp)
    subprocess.run(bash_command0, shell=True)
    subprocess.run(bash_command2, shell=True)
    with open("iter_file.txt", "w") as f:
        f.write (f"{1}")


