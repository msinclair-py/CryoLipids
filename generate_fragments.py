#!/usr/bin/env python
import json
from utilities import rtfParser
from molecular_graph import MolecularGraph

# OUTSTANDING ISSUE: since we start from fragments of size 3 there are many
#   fragments missing from any auxiliary paths including the sn2 tail. This
#   is likely trivial but should still be addressed in future update

# support for PI and Cardiolipin is pending as these are nontrivial
# to model even with our novel approach
lipids = ['POPE', 'POPC', 'POPG', 'POPS'] #, 'POPI24', 'PVCL2']

# get rtf for edge graph
parse = rtfParser()
rtfs = parse.get_rtfs

for lipid in lipids:
    print('\n---------------------')
    print(f'| Fragmenting {lipid}! |')
    print('---------------------')

    # generate adjacency lists
    graph = MolecularGraph(rtfs, lipid)
    G = graph.connectivity_graphs

    # using depth-first search, identify longest chain in molecule
    dfs = graph.longest_path
    print(f'\nLongest path identified:\n{"->".join(dfs)}')
    
    # from longest path, we can obtain a full library of fragments
    # of size [3, N] inclusive. this is nontrivial and is in fact
    # O(N*k^2) but exploiting lipid topology allows us to approach
    # this thoughtfully and obtain a fragment library in a fraction
    # of the time
    fragments = graph.fragment_lipid(dfs)
    print('\nFragments found:')
    for key, vals in fragments.items():
        print(f'Size {key}: {len(vals)}')

    json.dump(fragments, open(f'fragment_library/{lipid}.json', 'w'))
