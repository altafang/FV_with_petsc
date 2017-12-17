"""Go through all steps necessary to test and validate result."""

import subprocess
from make_inputs import make_sigma, make_source
from validate_results import validate_results

def run_test(variation='X'):
    
    if not ((variation == 'X') or (variation == 'Y')):
        print "Invalid dimension of variation, must be X or Y"
        return
        
    # Modify input.txt file
    subprocess.call("cp input_variation" + variation + ".txt input.txt", shell=True)

    command = "cp ../../bin/solve_poisson_2D ."
    subprocess.call(command, shell=True)
    print command

    # Make inputs
    print "Generating input files..."
    make_sigma(plot=False)
    make_source()

    # Run solve
    command = "mpiexec -np 1 ./solve_poisson_2D"
    print command
    print "Solving..."
    subprocess.call(command, shell=True)

    # Validate
    print "Validating..."
    validate_results(plot=False, variation=variation)
    
if __name__ == '__main__':
    run_test(variation='X')
    run_test(variation='Y')