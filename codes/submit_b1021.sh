#!/bin/bash
#SBATCH -A b1021                # Allocation
#SBATCH -p buyin                # Queue
#SBATCH -t 48:00:00             # Walltime/duration of the job
#SBATCH -N  1                   # Number of Nodes
#SBATCH --ntasks-per-node=28    # Number of Cores (Processors)
#SBATCH --job-name="LL_C16KK_NaCl10"    # Name of job

source /projects/b1021/Felipe/programs/bin/GMXRC

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt -ntmpi 4 -ntomp 7
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -n index.ndx -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -ntmpi 4 -ntomp 7
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -n index.ndx -p topol.top -o md.tpr
gmx mdrun -deffnm md -ntmpi 4 -ntomp 7
