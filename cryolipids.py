#!/usr/bin/env python
import MDAnalysis as mda
from MDAnalysis.analysis import align
from modeling import Lipid, Template
from molecular_graph import MolecularGraph
import numpy as np
from utilities import PDB, rtfParser
import tomllib

config = tomllib.load(open('config.toml','rb'))

# DELETE THIS LATER
structure = 'toy_models/model1.pdb'
#

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
pdb.write_to_pdb_file(pdb.contents)

# read in rtf files
parse = rtfParser()
rtfs = parse.get_rtfs

# NOTE: this will need to be reworked in order to be done
#   on the lipids from the above PDB() class structure

# read in lipid(s) to be modeled
molecule = Lipid(structure, 
                5, 
                rtfs['POPC'], 
                current_restype='POV')

# identify terminal atoms
terminal_atoms = molecule.get_terminal_atoms()

# get coords for stapling
tail_map = {'sn1': 'C2', 'sn2': 'C3'}
for (tail, terminus) in terminal_atoms.items():
        cur = f'{tail_map[tail]}{list(terminus.keys())[0]}'
        prev = f'{tail_map[tail]}{list(terminus.keys())[0] - 1}'
        prev_coords = molecule.get_coord(prev)
        vector_ref = np.array([list(terminus.values())[0], prev_coords])
        rtf_lipid = Template(f'lipids_from_rtf/{config["system"]["lipid"]}.pdb', 
                             config['system']['lipid'])
        vector_comp = rtf_lipid.atomic_coordinates([cur, prev])
        new_tail_names, new_tail_coords = rtf_lipid.missing_atoms(cur)
        new_tail_coords = molecule.staple_tail(vector_ref, vector_comp, new_tail_coords)
        new_tail_vec = [new_tail_coords[-1, :], new_tail_coords[0, :]]
        z = 1 # how do we decide 1 vs -1??
        new_tail_coords = molecule.staple_tail(np.array([[0,0,0], [0,0,z]]), 
                                               new_tail_vec, new_tail_coords)
        
        molecule.add_to_pdb(new_tail_names, new_tail_coords)
        molecule.write_to_pdb_file(molecule.pdb_contents)
        print(molecule.pdb_contents)

modeled_atoms = molecule.extract_coordinates()
#print(modeled_atoms)

# do minimization?

# check for protein conflict

# fix conflict

# do minimization?
