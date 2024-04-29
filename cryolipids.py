#!/usr/bin/env python
from collision_detection import CollisionDetector
from minimizer import VacuumSimulator, ImplicitSolventSimulator
from modeling import Lipid
from utilities import PDB, rtfParser
import tomllib

config = tomllib.load(open('config.toml','rb'))

# DELETE THIS LATER
structure = 'toy_models/model1.pdb'
#

# read in cryo em model file and process
pdb = PDB('misc/bmrcd.pdb', [5])
pdb.write_to_pdb_file(pdb.contents)

# load rtf files
rtfs = rtfParser()

# read in lipid(s) to be modeled
lipid = Lipid(structure, 
              5, 
              rtfs['POPC'], 
              current_restype='POV')

# model lipid
lipid.model()

# check for protein conflict; octree?
collision_detector = CollisionDetector(protein, lipid, method=0)
collision_detector.query_points()

# fix conflict


# do vacuum minimization
minimizer = VacuumSimulator(f'processed_{name}.pdb')

# do implicit solvent minimization and relaxation
minimizer = ImplicitSolventSimulator('vacuum_min.pdb')