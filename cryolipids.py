#!/usr/bin/env python
from modeling import Lipid
from utilities import PDB
import tomllib

config = tomllib.load(open('config.toml','rb'))

# DELETE THIS LATER
structure = 'toy_models/model1.pdb'
#

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
pdb.write_to_pdb_file(pdb.contents)

# read in lipid(s) to be modeled
molecule = Lipid(structure, 
                5, 
                rtfs['POPC'], 
                current_restype='POV')

molecule.model()

# do minimization?


# check for protein conflict; octree?


# fix conflict


# do minimization?
