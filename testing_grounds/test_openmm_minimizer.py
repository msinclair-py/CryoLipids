import openmm as mm
import openmm.app as app
from openmm.unit import *
from openmm.app import *
from openmm import *
import os
from sys import stdout, exit, stderr

print('Loading CHARMM PDB files ...')
# inpcrd = AmberInpcrdFile('/home/ss171/micromamba/envs/cryo/share/openmm/examples/input.inpcrd')
# prmtop = AmberPrmtopFile('/home/ss171/micromamba/envs/cryo/share/openmm/examples/input.prmtop', periodicBoxVectors=inpcrd.boxVectors)
charmm_input_pdb = 'test_lipid.pdb'
amber_output_pdb = 'test_lipid_amber.pdb'
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
command = f'amber4pdb -i {charmm_input_pdb} -o {amber_output_pdb}'

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
# system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# # Trying to setup system from AMBER formatted pdb
# pdb = PDBFile('test_lipid_ambpdb.pdb')
# forcefield = ForceField('amber14-all.xml')
# system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
#         nonbondedCutoff=1*nanometer, constraints=HBonds, implicitSolvent=OBC2)

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

'''
### Attempting to implement for charmm files

print('Loading CHARMM files...')
# params = CharmmParameterSet('toppar/par_all36_cgenff.prm', 
#                             'toppar/par_all36_lipid.prm', 
#                             'toppar/toppar_all36_lipid_cholesterol.prm',ls topp
#                             'toppar/toppar_water_ions.str')
toppar_dir = 'toppar'                              
param_files = ['par_all36m_prot.prm',
               'par_all36_cgenff.prm',
               'par_all36_lipid.prm',
               'toppar_all36_lipid_cholesterol_2.str',
               'toppar_water_ions.str']                                         

param_list = ['%s/%s'%(toppar_dir, param_file) for param_file in param_files]   
print(param_list)
params = CharmmParameterSet(*param_list)
print('Creating CHARMM psf and pdb...')
psf = CharmmPsfFile('step5_input.psf')
pdb = PDBFile('step5_input.pdb')
output_prefix = 'mineq' 
# make sure n_run_steps is divisible by chk_freq                                
chk_freq    = 250000                                                            
dcd_freq    = 5000                                                              
state_freq  = 5000    
  
psf.setBox(64.91400146484375*nanometer, 64.73699951171875*nanometer, 85.90900039672852*nanometer)

print('Creating CHARMM system...')
system = psf.createSystem(params, \
                          nonbondedMethod=PME,
                          nonbondedCutoff=1.2*nanometer,
                          switchDistance=1.0*nanometer,
                          constraints=HBonds,
                          rigidWater=True)

print('Starting engine...')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('%s.dcd'%(output_prefix), dcd_freq, enforcePeriodic
    Box=True))
# simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)

'''