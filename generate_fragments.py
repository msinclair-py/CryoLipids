#!/usr/bin/env python
import persim
import ripser
from file_utilities import PDB, rtfParser

min_frag = 5 # minimum number of atoms in a fragment for consideration
             # this number should correlate with the smallest collection of
             # modelled atoms we can have any confidence in completing
lipids = ['POPE', 'POPC', 'POPG', 'POPS'] #, 'POPI24', 'PVCL2']

# get rtf for edge graph
parse = rtfParser()
rtfs = parse.get_rtfs
graphs = parse.connectivity_graphs

for lipid in lipids:
    # get coords from lipids/*pdb
    lip = PDB(f'fragment_library/lipids/{lipid}', resname=lipid)
    coords = lip.contents

    # generate different length fragments
    heavy_atoms = [x for x in list(set().union(*graphs[lipid])) if 'H' not in x]
    connections = [c for c in graphs[lipid] if 'H' not in ''.join(c)]

    fragments = {2: connections}
    for n in range(3, len(heavy_atoms) - 1):
        fragments[n] = list()
        for frag in fragments[n-1]:
            for bond in connections:
                if frag.intersection(bond):
                    fragments[n].append(frag | bond)

    print(fragments)



    # iterate through all atoms to generate possible paths
    #for atom in all_atoms:
    #    conn = [connect for connect in connections if atom in connect]

    # write out to different fragment directories

    # PH diagrams

    # wasserstein distances
