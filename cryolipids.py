#!/usr/bin/env python
import argparse
import json
from utilities import PDB, rtfParser
from modeling import Lipid
from molecular_graph import MolecularGraph

parser = argparse.ArgumentParser(description = '')

parser.add_argument('structure', help='PDB to be modeled, must include \
        filepath if not in this directory. Output will be in same \
        directory as input.')
parser.add_argument('resids', help='ResID or ResIDs to be modeled. Note \
        that these must all be modeled the as the same lipid type at \
        this time. Must be a string (in quotes) and space-delimited.')
parser.add_argument('lipid', help='Type of lipid to build in place of \
        Cryo-EM modeled lipid.', 
        choices=['POPC', 'POPE', 'POPS', 'POPA', 'POPG'])

args = parser.parse_args()

structure = args.structure
resids = args.resids
lipid = args.lipid

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
print(pdb.contents)
pdb.write_to_pdb_file(pdb.contents)

# read in rtf files
parse = rtfParser()
rtfs = parse.get_rtfs

# read in lipid(s) to be modeled
molecule = Lipid('toy_models/model1.pdb', 
                5, 
                rtfs['POPC'], 
                graph._edges['POPC'], 
                current_restype='POV')

# generate persistence diagrams for each `lipid` fragment of size `k`
rtf_lipid = 
k = 

fragments = json.load('fragment_library/{lipid}.json')
diagrams = dict()
for frags in fragments[k]:
    for frag in frags:
        coords = rtf_lipid.atomic_coordinates(frag)
        diagram = graph.get_diagram(coords)
        print(diagram)

# identify fragment


# model

