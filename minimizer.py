import openmm as mm
import openmm.app as app
from openmm.unit import *
import os

class VacuumSimulator:
    """
    Simulation object which performs vacuum minimization. Can optionally tune any
    simulation parameters; however, the default settings should be appropriate for
    nearly any system.
    """
    def __init__(self, structure, output=os.getcwd(),
                 forcefield='amber14/protein.ff14SB.xml', 
                 temp=300 * kelvin, press=1 * bar, nonbondedMethod=app.NoCutoff, 
                 constraints=app.HBonds, collision_freq=1 / picosecond,
                 timestep=0.002 * picosecond, platform='CUDA'):
        self.structure = structure
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
        
    def set_integrator(self):
        pdb = app.PDBFile(self.structure)
        if self.solvent is not None:
            forcefield = app.ForceField(self.forcefield, self.solvent)
        else:
            forcefield = app.ForceField(self.forcefield)
            
        self.system = forcefield.createSystem(pdb.topology, nonbondedMethod=self.nonbondedMethod,
                                              constraints=self.constraints)
        self.integrator = mm.LangevinMiddleIntegrator(self.temp, self.collision_freq, self.ts)
        self.topology = pdb.topology
        self.init_pos = pdb.positions
        
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