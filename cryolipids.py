#!/usr/bin/env python
from utilities import PDB, rtfParser, MolecularGraph
from modeling import Lipid

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
pdb.write_to_pdb_file(pdb.contents)

# read in rtf files
parse = rtfParser()
rtfs = parse.get_rtfs

# generate molecular graph object
graph = MolecularGraph(rtfs)
graph_dict = graph.graphs
print(graph_dict['POPE'].edges)
#print(graph.adjacency())
#graph.connectivity_graphs
#adj_list = graph.adjacency_list

# read in lipid(s) to be modeled
#lipid = Lipid('toy_models/model1.pdb', 5, rtfs['POPC'], graph._edges['POPC'], restype='POV')
#lipid.model()
