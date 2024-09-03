#!/bin/bash

VMD="/Applications/VMD\\ 1.9.4a57-arm64-Rev12.app/Contents/vmd/vmd_MACOSXARM64"
LIPIDS=("POPE" "POPC" "POPG" "POPS" "POPI24" "POPA")

for lipid in ${LIPIDS[@]}; do
    cp phos.pdb temp.pdb
    sed -i "s/RES_A/$lipid/" temp.pdb
    $VMD -dispdev text -e gen_lipid.tcl -args $lipid
done
