import openmm as mm
import openmm.app as app
from openmm.unit import *
from openmm.app import *
from openmm import *
import os, shutil, tempfile
from sys import stdout, exit, stderr
import subprocess

class Simulator:
    """
    Simulation object which performs vacuum minimization. Can optionally tune any
    simulation parameters; however, the default settings should be appropriate for
    nearly any system.

    Because AMBER has limited lipid support, users will need to ensure that any lipid
    they wish to model have existing parameter information. This can be found at:
    https://ambermd.org/tutorials/advanced/tutorial16/#Lipid14
    
    To convert the charmm input files to amber input files, we will use a series of
    packages that come with ambertools. It's nice that we are just working with lipids
    because that means, there are no strange patches to consider when changing atom types
    and names.
    
    We will use, in order, the following ambertools packages:
    a) charmmlipid2amber.py -i input_structure.pdb [-c substitution_definitions.csv] -o output_structure.pdb
    b) pdb4amber -i {input file} -o {output file}
    c) tleap -s -f {tleap.in file} > {tleap.out log file}
    
    (A) and (B) should be fairly straight forward and do not require much tinkering,
    however, (C) will require some specifics to be included that were not initially 
    obvious to me.

    The only nonbonded methods that are supported with implicit solvent are NoCutoff (the default), 
    CutoffNonPeriodic , and CutoffPeriodic. If you choose to use a nonbonded cutoff with 
    implicit solvent, it is usually best to set the cutoff distance larger than is typical 
    with explicit solvent.
    """
    def __init__(self, structure, output=os.getcwd(),
                 forcefield=['amber14-all.xml', 'amber14/protein.ff14SB.xml'], 
                 temp=300 * kelvin, press=1 * bar, nonbondedMethod=app.NoCutoff, 
                 constraints=app.HBonds, collision_freq=1 / picosecond,
                 timestep=0.002 * picosecond, platform='CUDA', solvent=None):
        self.charmm_structure = structure
        self.output = output
        self.forcefield = forcefield
        self.temp = temp
        self.pressure = press
        self.nonbondedMethod = nonbondedMethod
        self.constraints = constraints
        self.collision_freq = collision_freq
        self.ts = timestep
        self.platform = platform # need to do a gpu check here
        self.solvent = solvent # allows us to inherit methods more cleanly
        self.rst7 = 'amber_file.rst7'
        self.prmtop = 'amber_file.prmtop'
        self.tleap_template = 'tleap.in'
        self.amber_tmp_file = 'amber_format.pdb'
        # self.output = '' # this is for all output in commands list
        # self.minimized = 'minimize_run.pdb'
        # self.vacuumed = 'vacuum_run.pdb'
        # self.implicited = 'implicit_run.pdb'
        self.simulation = None # use this to access simulation object
    
    def prep(self):
        # Convert CHARMM lipid naming to AMBER convention
        commands = [f"charmmlipid2amber.py -i {self.charmm_structure} -o renamed_lipids.pdb -c charmmlipid2.amber.csv",
                    f"pdb4amber -i renamed_lipids.pdb -o {self.amber_tmp_file} -y",
                    f'sed -i "s/CD  ILE/CD1 ILE/" {self.amber_tmp_file}']
        
        # Convert CHARMM PDB file to AMBER formatting
        subprocess.run(commands[0], shell=True, capture_output=True, text=True)
        subprocess.run(commands[1], shell=True, capture_output=True, text=True)
        
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.tleap') as temp_file:
            self.tleap_conf = temp_file.name
            
            # Read the original tleap file and replace placeholders
            with open(self.tleap_template, 'r') as original_file:
                content = original_file.read()
                content = content.replace('{amber_format}', self.amber_tmp_file)
                content = content.replace('{prmtop_file}', self.prmtop)
                content = content.replace('{rst7_file}', self.rst7)
            
            # Write the modified content to the temporary file
            temp_file.write(content)
            temp_file.flush()  # Ensure content is written to disk
        
            # Use the modified temporary file with tleap within the block
            try:
                command = f'tleap -s -f {self.tleap_conf} > tleap.log'
                subprocess.run(command, shell=True, capture_output=True, text=True)
            except Exception as e:
                print(f'Error running tleap: {e}')
        
        # The temporary file is automatically deleted after the with block ends
        print('\n Prep done :) \n')
        
    def minimize(self, solvent=None):

        self.solvent = solvent 
        
        inpcrd = AmberInpcrdFile(self.rst7)
        prmtop = AmberPrmtopFile(self.prmtop)
        
        if self.solvent is None:
            system = prmtop.createSystem(constraints=HBonds)

            integrator = LangevinMiddleIntegrator(self.temp, self.collision_freq, self.ts)
            self.simulation = Simulation(prmtop.topology, system, integrator)
            self.simulation.context.setPositions(inpcrd.positions)
            self.simulation.minimizeEnergy()
            
            with open('vacuum_minimized.pdb', 'w') as output:
                state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                PDBFile.writeFile(self.simulation.topology, state.getPositions(), output)
                self.simulation.saveState('vacuum_minimized.xml')

            print('\n Vacuum minimizer done :) \n')

        elif solvent == 'implicit':
            # forcefield = ForceField(*self.forcefield, self.solvent)
            system = prmtop.createSystem(constraints=HBonds, implicitSolvent=GBn2) # implicit solvent

            integrator = LangevinMiddleIntegrator(self.temp, self.collision_freq, self.ts)
            self.simulation = Simulation(prmtop.topology, system, integrator)
            self.simulation.loadState('vacuum_minimized.xml')
            self.simulation.context.setPositions(inpcrd.positions)
            self.simulation.minimizeEnergy()

            with open('implicit_minimized.pdb', 'w') as output:
                state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                PDBFile.writeFile(self.simulation.topology, state.getPositions(), output)
                self.simulation.saveState('implicit_minimized.xml')
            print('\n Implicit minimizer done :) \n')


    def vacuum_run(self):
        self.simulation.reporters.append(DCDReporter('vacuum.dcd', 1000))
        self.simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        self.simulation.step(10000)
        
        # Test write PDB to ensure that the minimize method is not killing my PDB coordinates
        # openmm github thread: https://github.com/openmm/openmm/issues/3635
        with open('vacuumed.pdb', 'w') as output:
            state = self.simulation.context.getState(getPositions=True, getEnergy=True)
            PDBFile.writeFile(self.simulation.topology, state.getPositions(), output)
            self.simulation.saveState('vacuum_run.xml')
            
    def implicit_run(self):
        self.simulation.reporters.append(DCDReporter('implicit.dcd', 1000))
        self.simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                potentialEnergy=True, temperature=True))
        self.simulation.step(10000)
        
        # Test write PDB to ensure that the minimize method is not killing my PDB coordinates
        # openmm github thread: https://github.com/openmm/openmm/issues/3635
        with open('implicited.pdb', 'w') as output:
            state = self.simulation.context.getState(getPositions=True, getEnergy=True)
            PDBFile.writeFile(self.simulation.topology, state.getPositions(), output)
            self.simulation.saveState('implicit_run.xml')

