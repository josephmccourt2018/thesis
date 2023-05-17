#!/bin/bash

basedir=$(pwd)
for state in 10 20 30 40 50 60 70 80 90 100

do
    cd $state
    mkdir combined
    cd combined
    python ../../area_per_lipid_trajectory.py
    gmx make_ndx -f ../npt.gro -o dens_groups.ndx<<EOF
r LYS | r LSN
r C16
q
EOF
    #gmx make_ndx -f npt.gro -o c_alpha.ndx<<EOF
#r LYS | r LSN & a CA 
#q
#EOF

    gmx make_ndx -f ../npt.gro -o order_index.ndx<<EOF
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
del 0
r C16 & a CA
r C16 & a CB
r C16 & a CC
r C16 & a CD
r C16 & a CE
r C16 & a CF
r C16 & a CG
r C16 & a CH
r C16 & a CI
r C16 & a CJ
r C16 & a CK
r C16 & a CL
r C16 & a CM
r C16 & a CN
r C16 & a CO
r C16 & a C
q
EOF

    echo "18" | gmx density -f ../combined.xtc -s ../md_semi.tpr -n dens_groups.ndx -o dens_head.xvg -sl 100
    echo "19" | gmx density -f ../combined.xtc -s ../md_semi.tpr -n dens_groups.ndx -o dens_tails.xvg -sl 100
    echo "13" | gmx density -f ../combined.xtc -s ../md_semi.tpr -n dens_groups.ndx -o dens_water.xvg -sl 100
    echo "16" | gmx density -f ../combined.xtc -s ../md_semi.tpr -n dens_groups.ndx -o dens_ions.xvg -sl 100
    gmx order -f ../combined.xtc -o order.xvg -s ../md_semi.tpr -n order_index.ndx
    #gmx rdf -f md_semi.xtc -n c_alpha.ndx -o rdf_chiral_carbon.xvg -ref 18 -sel 18
    #gmx energy -f md_semi.edr -o LJ_sr_pp.xvg<<EOF
#46
#0
#EOF
    #gmx energy -f md_semi.edr -o coul_sr_pp.xvg<<EOF
#45
#0
#EOF
    
    cd $basedir
done
