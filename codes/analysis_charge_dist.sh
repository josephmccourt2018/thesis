#!/bin/bash

basedir=$(pwd)
for state in 10 20 30 40 50 60 70 80 90 100

do
    cd 10nm_box
    cd $state
    rm *\#
    gmx make_ndx -f md_semi_cont.gro -o lys_groups_cont.ndx<<EOF
r LYS & a CA
r LYS & a NZ
r LSN & a NZ
r LSN | r LYS & a NZ
q
EOF

    gmx rdf -f md_semi_cont.xtc -n lys_groups_cont.ndx -o rdf_lys_chiral_carbon_cont.xvg -ref 18 -sel 18
    gmx rdf -f md_semi_cont.xtc -n lys_groups_cont.ndx -o rdf_lys_Nz_cont.xvg -ref 19 -sel 19
    gmx rdf -f md_semi_cont.xtc -n lys_groups_cont.ndx -o rdf_lsn_Nz_cont.xvg -ref 20 -sel 20
    gmx rdf -f md_semi_cont.xtc -n lys_groups_cont.ndx -o rdf_lys_lsn_Nz_cont.xvg -ref 21 -sel 21
    
    cd $basedir
done
