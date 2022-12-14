#!/usr/bin/env python
import json
from utilities import rtfParser
from molecular_graph import MolecularGraph

lipids = ['POPE'] #, 'POPC', 'POPG', 'POPS'] #, 'POPI24', 'PVCL2']

# get rtf for edge graph
parse = rtfParser()
rtfs = parse.get_rtfs

for lipid in lipids:
    print(f'Fragmenting {lipid}!')
    # generate adjacency lists
    graph = MolecularGraph(rtfs, lipid)
    G = graph.connectivity_graphs
    
    dfs = graph.longest_path
    print(f'Longest path identified:\n{"->".join(dfs)}')
    
    fragments = graph.fragment_lipid(dfs)
    #print(fragments[15])


    #fragments = graph.fragment_lipid
    #json.dump(fragments, open(f'fragment_library/{lipid}.json', 'w'))
