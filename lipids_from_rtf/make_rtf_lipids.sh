#!/bin/bash

LIPIDS=("POPE" "POPC" "POPG" "POPS" "POPI24")

for lipid in ${LIPIDS[@]}; do
    cp phos.pdb temp.pdb
    sed -i "s/RES_/$lipid/" temp.pdb
    /Common/linux/bin/vmd -dispdev text -e gen_lipid.tcl -args $lipid
done
