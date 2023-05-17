#!/bin/bash
#SBATCH -A b1030                # Allocation
#SBATCH -p buyin                # Queue
#SBATCH -t 48:00:00             # Walltime/duration of the job
#SBATCH -N  1                   # Number of Nodes
#SBATCH --ntasks-per-node=28    # Number of Cores (Processors)
#SBATCH --gres=gpu:p100:2 
#SBATCH --job-name="C16K_bilayer"    # Name of job

source /projects/b1021/Felipe/programs/bin/GMXRC

gmx grompp -f npt.mdp -c minim.gro -r minim.gro -t minim.cpt -n index.ndx -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -ntmpi 4 -ntomp 7 -gpu_id 0011 -maxh 48
gmx grompp -f prod_semi.mdp -c npt.gro -t npt.cpt -n index.ndx -p topol.top -o md.tpr
gmx mdrun -deffnm md -ntmpi 4 -ntomp 7 -maxh 36 -gpu_id 0011 -maxh 48
