#!/usr/bin/env python
from collision_detection import Repairer
from minimizer import Simulator
from modeling import Lipid
from utilities import PDB, unpack_lipids
import tomllib

config = tomllib.load(open('config.toml','rb'))
input_name = config['system']['name']
output_name = config['system']['output']

lipid_resids, lipid_types, modeled_resnames = unpack_lipids(config)

# read in cryo em model file and process
print('Starting setup for modeling lipids')
pdb = PDB(input_name, output_name, lipid_resids, modeled_resnames)
incomplete_lipids = pdb.write_to_pdb_file(pdb.contents)
protein_coords = [[float(i) for i in line[31:54].strip().split()] 
                  for line in pdb.protein]

new_coords = dict()

# modeling individual lipids
print('Modeling lipids individually')
for id, lip in config['lipids'].items():
    print(id)
    # read in lipid(s) to be modeled
    lipid = Lipid(incomplete_lipids, **lip) 
    print('test 1')
    # model lipid
    lipid.model()
    print('test 2')
    # check for atomic clashes and repair accordingly
    repair = Repairer(lipid, protein_coords, grid_spacing=1.5)
    print('repaired')
    repair.check_collisions()
    print('no collisions')
    new_coords[id] = repair.get_new_coords()

# output final static model
pdb.merge_final_pdb(new_coords)

shashank_test = 'shashank_lipid_test.pdb'
if config['minimize']['vacuum']:
    # do vacuum minimization
    # minimizer = VacuumSimulator(f'{output_name}.pdb')
    print(f'working on file: {shashank_test}')
    minimizer = Simulator(f'{shashank_test}')
    minimizer.prep()
    minimizer.minimize()
    
    if config['minimize']['implicit_solvent']:
        # do implicit solvent minimization and relaxation
        # minimizer = ImplicitSolventSimulator('vacuum_min.pdb')
        minimizer.minimize(solvent='implicit')
