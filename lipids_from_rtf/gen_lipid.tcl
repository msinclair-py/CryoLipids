package require psfgen

set RES [lindex $argv 0]

topology /Scr/msincla01/c36/top_all36_lipid.rtf
topology /Scr/msincla01/c36/toppar_all36_lipid_inositol.str

segment LIP {
    pdb temp.pdb
}

coordpdb temp.pdb LIP
guesscoord
regenerate angles dihedrals

writepdb $RES.pdb

exit
