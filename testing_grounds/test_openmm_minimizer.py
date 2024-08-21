import openmm as mm
import openmm.app as app
from openmm.unit import *
from openmm.app import *
# from openmm import *
import os
from sys import stdout, exit, stderr

print('Loading CHARMM PDB files ...')
charmm_input_pdb = 'test_lipid.pdb'
amber_input_pdb = 'test_lipid_amber.pdb'

print('Converting to AMBER PDB files ...')

'''
To convert the charmm input files to amber input files, we will use a series of
packages that come with ambertools. It's nice that we are just working with lipids
because that means, there are no strange patches to consider when changing atom types
and names.

We will use, in order, the following ambertools packages:
a) charmmlipid2amber.py
b) pdb4amber -i {input file} -o {output file}
c) tleap -s -f {tleap.in file} > {tleap.out log file}

(A) and (B) should be fairly straight forward and do not require much tinkering,
however, (C) will require some specifics to be included that were not initially 
obvious to me.
- 
'''
# Construct the command using f-strings
command = f'pdb4amber -i {charmm_input_pdb} -o {amber_output_pdb}'

# Run the command
# result = subprocess.run(command, shell=True, capture_output=True, text=True)

# Trying to setup system using prmtop/rst7 files
# do not try to generate prmtop/inpcrd files as they don't seem
# to provide the information necessary for openmm to start

inpcrd = AmberInpcrdFile('test_lipid.rst7')
prmtop = AmberPrmtopFile('test_lipid.prmtop')
forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml', 'implicit/gbn2.xml')

print('Building system ...')
system = prmtop.createSystem(constraints=HBonds, implicitSolvent=GBn2)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()


# Test write PDB to ensure that the minimize method is not killing my PDB coordinates
# openmm github thread: https://github.com/openmm/openmm/issues/3635
with open('minimized.pdb', 'w') as output:
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    PDBFile.writeFile(simulation.topology, state.getPositions(), output)

simulation.reporters.append(DCDReporter('output.dcd', 1))
# simulation.reporters.append(PDBReporter('output.pdb', 5))
simulation.reporters.append(StateDataReporter(stdout, 5, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(500)
