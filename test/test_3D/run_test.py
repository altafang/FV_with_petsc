"""Go through all steps necessary to test and validate result."""

import subprocess
from make_inputs import make_sigma, make_source
from validate_results import validate_results

command = "cp ../../bin/solve_poisson_3D ."
subprocess.call(command, shell=True)
print command

# Make inputs
print "Generating input files..."
make_sigma(plot=False)
make_source()

# Run solve
command = "mpiexec -np 1 ./solve_poisson_3D"
print command
print "Solving..."
subprocess.call(command, shell=True)

# Validate
print "Validating..."
validate_results(plot=False)