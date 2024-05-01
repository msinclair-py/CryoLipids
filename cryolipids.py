#!/usr/bin/env python
from collision_detection import CollisionDetector
from minimizer import VacuumSimulator, ImplicitSolventSimulator
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
lipid = Lipid(structure, 
              5, 'POPC',
              current_restype='POV')

# model lipid
lipid.model()

# check for protein conflict; octree?
protein_coords = [[float(i) for i in line[32:54].strip().split()] 
                  for line in pdb.protein]

collision_detector = CollisionDetector(protein_coords, 
                                       lipid.extract_coordinates(), 
                                       method=0)

collisions = collision_detector.query_points()
print(collisions)

# fix conflict
if any(collisions):
    pass

# do vacuum minimization
minimizer = VacuumSimulator(f'processed_{name}.pdb')

# do implicit solvent minimization and relaxation
minimizer = ImplicitSolventSimulator('vacuum_min.pdb')