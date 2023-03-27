#!/usr/bin/env python
import json
from modeling import Template
from molecular_graph import PersistentHomology
from optimization import Optimizer

lipids = ['POPE', 'POPC', 'POPG', 'POPS']
k_upper_bounds = [48, 51, 50, 51] # N_{heavy atoms} - 1
sample_size = 5 # no. fragments for each size `k` to sample

# generate persistence diagrams for each `lipid` fragment of size `k`
for (lipid, k_upper) in zip(lipids, ks):
    fragments = json.load(open(f'fragment_library/{lipid}.json', 'r'))

    diagrams = dict()
    for i, frag in enumerate(fragments[k]):
        coords = rtf_lipid.atomic_coordinates(frag)
        diagram = PersistentHomology.get_diagram(coords)
    
        diagrams[i] = diagram

    # sample fragments
    test_dict = dict()
    for k in range(3, k_upper):
        test_dict[k] = {}
        test_idxs = np.random.choice(len(fragments[k]), 
                                    size=sample_size, 
                                    replace=False)
        temp = {idx: frag[idx] for idx, frag in enumerate(fragments[k])}

        # introduce noise to data
        for i, fragment in enumerate(fragments[k][test_idxs]):
            fragment = add_noise(fragment)
            test_dict[k].update({i: fragment})

    
    # validate linear combination
    
    # identify fragment
    PH = PersistentHomology(modeled_atoms, diagrams)
    fragment_index = PH.best_fit
    atom_names = fragments[k][fragment_index]
