#!/bin/bash
#SBATCH --job-name=p2_eps01_no_skip
#SBATCH --output=p2_eps01_no_skip_%j.out
#SBATCH --error=p2_eps01_no_skip_%j.err
#SBATCH --partition=long
#SBATCH --nodes=6                                                    
#SBATCH --ntasks-per-node=24 
#SBATCH --ntasks=144                                  
#SBATCH --exclusive
#SBATCH --time=3-00:00:00                             


source /etc/profile.d/modules.sh
module purge
module load cmake gcc openmpi

#mpirun -n 144 ./poet dolo_fgcs_3.R dolo_fgcs_3.qs2 dolo_only_pqc
mpirun -n 144 ./poet --interp --rds dolo_fgcs_3_rt.R dolo_fgcs_3.qs2 p2_eps01_no_skip
#mpirun -n 144 ./poet --interp  barite_fgcs_4_new/barite_fgcs_4_new_rt.R barite_fgcs_4_new/barite_fgcs_4_new.qs2 barite