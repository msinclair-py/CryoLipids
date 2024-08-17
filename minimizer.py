import openmm as mm
import openmm.app as app
from openmm.unit import *
from openmm.app import *
from openmm import *
import os
from sys import stdout, exit, stderr

class VacuumSimulator:
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
                 timestep=0.002 * picosecond, platform='CUDA'):
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
        self.solvent = None # allows us to inherit methods more cleanly
        self.tleap_conf = 
        self.tleap_log =
        self.rst7 = ''
        self.prmtop = ''

        # Convert CHARMM lipid naming to AMBER convention
        command = f'charmmlipid2amber.py -i {self.charmm_structure} \
                    -o {output}/renamed_lipids.pdb'
        try:
            subprocess.run(command, shell=True, capture_output=True, text=True)
        except:
            print('Fixme -- write a specific exception')

        # Convert CHARMM PDB file to AMBER formatting
        command = f'pdb4amber -i {output}/renamed_lipids.pdb -o {output}/amber_format.pdb'
        try:
            subprocess.run(command, shell=True, capture_output=True, text=True)
        except:
            print('Fixme -- write a specific exception')

        # Generate inpcrd and rst7 files with tleap
        # files are named amber_lipids.{inpcrd,prmtop}
        command = f'tleap -s -f {self.tleap_conf} > {self.tleap_log}'

    def set_integrator(self):
        rst7 = AmberInpcrdFile(self.rst7)
        prmtop = AmberPrmtopFile(self.prmtop)
        if self.solvent is not None:
            forcefield = ForceField(*self.forcefield, self.solvent)
        else:
            forcefield = ForceField(self.forcefield, 'implicit/gbn2.xml')

        # nonbonded method 
        # self.system = forcefield.createSystem(pdb.topology, nonbondedMethod=self.nonbondedMethod,
        #                                       constraints=self.constraints)
        self.system = prmtop.createSystem(constraints=self.constraints, implicitSolvent=GBn2)
        self.integrator = mm.LangevinMiddleIntegrator(self.temp, self.collision_freq, self.ts)
        self.topology = prmtop.topology
        self.init_pos = inpcrd.positions
        
    def minimize(self, positions):
        self.set_integrator()
        mm.Platform.getPlatformByName(self.platform)
        simulation = app.Simulation(self.topology, self.system, self.integrator)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        
        self.write_to_pdb(simulation.topology, 
                          simulation.context.getState(getPositions=True).getPositions(),
                          os.path.join(self.output, 'min.pdb'))
    
    def propagate_dynamics(self, positions, n_steps, out_name, save_rate, velocities=None):
        self.set_integrator()
        mm.Platform.getPlatformByName(self.platform)
        simulation = app.Simulation(self.topology, self.system, self.integrator)
        simulation.context.setPositions(positions)
        if velocities:
            simulation.setVelocities(velocities)
        simulation.reporters.append(app.DCDReporter(out_name, save_rate))
        simulation.step(n_steps)
        
        self.write_to_pdb(simulation.topology,
                          simulation.context.getState(getPositions=True).getPositions(),
                          os.path.join(self.output, out_name))
        
    @staticmethod
    def write_to_pdb(top, positions, out):
        with open(out, 'w') as outfile:
            app.pdbfile.PDBFile.writeFile(top, positions, outfile)
      
        
class ImplicitSolventSimulator(VacuumSimulator):
    """
    Simulation object for implicit solvent simulation. Defaults to GBN2
    implicit solvent. For other options see OpenMM documentation:
    """
    def __init__(self, structure, output=os.getcwd(),
                 forcefield='amber14/protein.ff14SB.xml', 
                 temp=300 * kelvin, press=1 * bar, nonbondedMethod=app.NoCutoff, 
                 constraints=app.HBonds, collision_freq=1 / picosecond,
                 timestep=0.002 * picosecond, platform='CUDA', 
                 solvent='implicit/gbn2.xml'):
        super().__init__(structure, output, forcefield, temp, press, nonbondedMethod,
                 constraints, collision_freq, timestep, platform)
        self.solvent = solvent