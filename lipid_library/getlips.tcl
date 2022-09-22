mol new membrane.psf
mol addfile membrane.pdb

for {set i 1} {$i < 19} {incr i} {
    set lip [atomselect top "segname MEMB and resid $i"]
    set name [lsort -u [$lip get resname]]
    set com [lindex [measure center [atomselect top "name P P1 P3 O3 and resid $i"]] 2]
    
    if {$com > 0} {        
        set leaf "upper"
    } else {
        set leaf "lower"
    }

    set cent [vecinvert [measure center $lip]]
    set movement [list [lindex $cent 0] [lindex $cent 1] 0]
    $lip moveby $movement
    $lip set resid 1

    $lip writepdb ${name}_${leaf}.pdb
    $lip writepsf ${name}.psf
}
