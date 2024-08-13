package require psfgen
resetpsf

set RES POPA
#set RES [lindex $argv 0]

topology ../rtf_files/top_all36_lipid.rtf
topology ../rtf_files/toppar_all36_lipid_inositol.str

segment LIP {
    pdb temp.pdb
}

coordpdb temp.pdb LIP
guesscoord
regenerate angles dihedrals

writepdb $RES.pdb

exit
