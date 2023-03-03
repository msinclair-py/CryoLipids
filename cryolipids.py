#!/usr/bin/env python
import argparse
import json
import MDAnalysis as mda
from MDAnalysis.analysis import align
from modeling import Lipid, Template
from molecular_graph import MolecularGraph, PersistentHomology
from utilities import PDB, rtfParser

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

# DELETE THIS LATER
structure = 'toy_models/model1.pdb'

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

modeled_atoms = molecule.extract_coordinates()

# generate persistence diagrams for each `lipid` fragment of size `k`
rtf_lipid = Template(f'lipids_from_rtf/{lipid}.pdb', lipid)
k = str(len(modeled_atoms))

fragments = json.load(open(f'fragment_library/{lipid}.json', 'r'))
diagrams = dict()
for i, frag in enumerate(fragments[k]):
    coords = rtf_lipid.atomic_coordinates(frag)
    diagram = PersistentHomology.get_diagram(coords)

    diagrams[i] = diagram

# identify fragment
PH = PersistentHomology(modeled_atoms, diagrams)
fragment_index = PH.best_fit
atom_names = fragments[k][fragment_index]

print(fragment_index, atom_names)

# align proper atoms to fragment
good_lipid = mda.Universe(rtf_lipid.filename)
good_CoM = good_lipid.select_atoms(atom_names).center_of_mass()
compare = good_lipid.select_atoms(atom_names).positions - good_CoM

ref = mda.Universe(structure)
ref_CoM = ref.select_atoms('all').center_of_mass()
reference = ref.select_atoms('all').positions - ref_CoM

rotation_matrix, _ = align.rotation_matrix(compare, reference)
good_lipid.atoms.translate(-good_CoM)
good_lipid.atoms.rotate(rotation_matrix)
good_lipid.atoms.translate(ref_CoM)

# check for protein clashes

# adjust structure accordingly

