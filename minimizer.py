import openmm as mm
import openmm.app as app
from openmm.unit import *

class VacuumSimulator:
    """
    
    """
    def __init__(self, structure, forcefield="amber14/protein.ff14SB.xml", 
                 temp=300 * kelvin, press=1 * bar, nonbondedMethod=app.NoCutoff, 
                 constraints=app.HBonds, collision_freq=1 / picosecond,
                 timestep=0.002 * picosecond, platform="CUDA"):
        self.structure = structure
        self.forcefield = forcefield
        self.temp = temp
        self.pressure = press
        self.nonbondedMethod = nonbondedMethod
        self.constraints = constraints
        self.collision_freq = collision_freq
        self.ts = timestep
        self.platform = platform # need to do a gpu check here
        
    def set_integrator(self):
        pdb = app.PDBFile(self.structure)
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
    
    def propagate_dynamics(self, positions, n_steps, out_name, save_rate, velocities=None):
        self.set_integrator()
        mm.Platform.getPlatformByName(self.platform)
        simulation = app.Simulation(self.topology, self.system, self.integrator)
        simulation.context.setPositions(positions)
        if velocities:
            simulation.setVelocities(velocities)
        simulation.reporters.append(app.DCDReporter(out_name, save_rate))
        simulation.step(n_steps)
      
        
class ImplicitSolventSimulator(VacuumSimulator):
    def __init__(self):
        super().__init__()
        pass