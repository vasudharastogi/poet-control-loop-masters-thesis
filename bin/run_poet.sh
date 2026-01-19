#!/bin/bash
#SBATCH --job-name=m_p1_50_1_3_1
#SBATCH --output=m_p1_50_1_3_1_%j.out
#SBATCH --error=m_p1_50_1_3_1_%j.err
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
mpirun -n 144 ./poet --interp dolo_fgcs_3_rt.R dolo_fgcs_3.qs2 m_p1_50_1_3_1
#mpirun -n 144 ./poet --interp  barite_fgcs_4_new/barite_fgcs_4_new_rt.R barite_fgcs_4_new/barite_fgcs_4_new.qs2 barite