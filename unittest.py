from minimizer import Simulator
import os, subprocess

# output_name = 'test_membrane.pdb'
output_name = 'test.pdb'

print('Creating simulator object ...', end='\n')
print('Original set of files in the system:', subprocess.run(['ls', '-ltr'], capture_output=True, text=True).stdout, end='\n')

minimizer = Simulator(output_name, unittest=True)

minimizer.prep()

minimizer.minimize()

minimizer.minimize(solvent='implicit')