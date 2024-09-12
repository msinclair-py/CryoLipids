#!/usr/bin/env python
from collision_detection import Repairer
from minimizer import Simulator
from modeling import Lipid
from utilities import PDB, unpack_lipids
import tomllib

config = tomllib.load(open('config.toml','rb'))
input_name = config['system']['name']
output_name = config['system']['output']

theta, min_attempts = config['repair'].values()

lipid_resids, lipid_types, modeled_resnames = unpack_lipids(config)

# read in cryo em model file and process
pdb = PDB(input_name, output_name, lipid_resids, modeled_resnames)
incomplete_lipids = pdb.write_to_pdb_file(pdb.contents)
protein_coords = [[float(i) for i in line[31:54].strip().split()] 
                  for line in pdb.protein]

new_coords = dict()
for id_, lip in config['lipids'].items():
    # read in lipid(s) to be modeled
    lipid = Lipid(incomplete_lipids, **lip) 

    # model lipid
    lipid.model()
    print('done modeling')

    # check for atomic clashes and repair accordingly
    repair = Repairer(lipid, protein_coords, theta, min_attempts, grid_spacing=1.5)
    repair.check_collisions()
    new_coords[id_] = repair.get_new_coords()
    print('done repairing')

# output final static model
pdb.merge_final_pdb(new_coords)

if config['minimize']['vacuum']:
    # do vacuum minimization
    minimizer = Simulator(output_name)
    minimizer.prep()
    minimizer.minimize()
    
    if config['minimize']['implicit_solvent']:
        # do implicit solvent minimization and relaxation
        minimizer.minimize(solvent='implicit')
