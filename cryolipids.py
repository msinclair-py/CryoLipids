#!/usr/bin/env python
from file_utilities import PDB, rtfParser
from modeling import Lipid

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
pdb.write_to_pdb_file(pdb.contents)

# read in rtf files
parse = rtfParser()
rtfs = parse.get_rtfs
graphs = parse.connectivity_graphs

# read in lipid(s) to be modeled
lipid = Lipid('toy_models/model1.pdb', 5, rtfs['POPC'], graphs['POPC'], restype='POV')
print(graphs)
lipid.model()
