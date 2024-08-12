#!/usr/bin/env python
from collision_detection import Repairer
#from minimizer import VacuumSimulator, ImplicitSolventSimulator
from modeling import Lipid
from utilities import PDB
import tomllib

config = tomllib.load(open('config.toml','rb'))

# DELETE THIS LATER
structure = 'toy_models/model1.pdb'
#

lipids_to_model = [5] #[5, 6, 7, 11] # list of resids

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', lipids_to_model)
pdb.write_to_pdb_file(pdb.contents) # writes initial lipid coords??

# read in lipid(s) to be modeled
lipid = Lipid(structure, 
              5, 'POPC',
              current_restype='POV')

# model lipid
lipid.model()

# check for protein conflicts
protein_coords = [[float(i) for i in line[31:54].strip().split()] 
                  for line in pdb.protein]

repair = Repairer(lipid, protein_coords, grid_spacing=1.5)
repair.check_collisions()
repair.write_pdb()

# do vacuum minimization
#minimizer = VacuumSimulator(f'processed_{name}.pdb')

# do implicit solvent minimization and relaxation
#minimizer = ImplicitSolventSimulator('vacuum_min.pdb')
