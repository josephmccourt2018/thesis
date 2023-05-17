#!/bin/bash

#basedir=$(pwd)
for state in 100
	     
do

    cd $state

    cd 4/bilayer
    cp ../../../prod_semi_long.mdp .
    gmx grompp -f npt_semi.mdp -c minim.gro -o npt_semi.tpr -n index.ndx -p topol.top
    gmx mdrun -deffnm npt_semi -ntmpi 4 -ntomp 7 -gpu_id 0011
    gmx grompp -f prod_semi.mdp -c npt_semi.gro -o md.tpr -n index.ndx -p topol.top
    gmx mdrun -deffnm md -ntmpi 4 -ntomp 7 -gpu_id 0011
    gmx grompp -f prod_semi_long.mdp -c md.gro -t md.cpt -n index.ndx -p topol.top -o md_2.tpr
    gmx mdrun -deffnm md_2 -ntmpi 4 -ntomp 7 -gpu_id 0011
    
    
    cd $basedir
    
done


#for the 0 charged
#cd ../4_with_unit_5/bilayer
#cp ../../charge_variation/prod_semi_long.mdp .
#gmx grompp -f prod_semi_long.mdp -c md.gro -t md.cpt -n index.ndx -p topol.top -o md_2.tpr
#gmx grompp -f prod_semi_long.mdp -c md.gro -n index.ndx -p topol.top -o md_2.tpr
#gmx mdrun -deffnm md_2 -ntmpi 4 -ntomp 7 -gpu_id 0011	  
