#!/usr/bin/env python
from utilities import PDB, rtfParser
from modeling import Lipid
from molecular_graph import MolecularGraph

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
print(pdb.contents)
pdb.write_to_pdb_file(pdb.contents)

# read in rtf files
parse = rtfParser()
rtfs = parse.get_rtfs

# generate molecular graph object
graph = MolecularGraph(rtfs)
graph.connectivity_graphs
adj_list = graph.adjacency_list

dummy = MolecularGraph.matt_algorithm(adj_list['POPE'], 6)
print(dummy)

# generate persistence diagrams for each fragment
rtf_lipid = Template(f'lipids_from_rtf/{lipid}', lipid)
diagrams = dict()
for k, fragment in fragments.items():
    for frag in fragment:
        coords = rtf_lipid.atomic_coordinates(frag)
        diagram = graph.get_diagram(coords)
        print(diagram)


# read in lipid(s) to be modeled
#lipid = Lipid('toy_models/model1.pdb', 5, rtfs['POPC'], graph._edges['POPC'], restype='POV')
#lipid.model()
