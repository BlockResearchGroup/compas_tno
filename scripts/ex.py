import os

from cvxpy import *
import sdpt3glue

# Declare variables:
x = Variable(1)

# Define objective:
obj = Minimize(2*x+10)

# Define constraints
constraints = [x >= 0]

# Construct the Cvxpy problem
problem = Problem(obj, constraints)
print(problem.solve())
print(x.value)

# Generate filenames
folder = '/Users/mricardo/Documents/MATLAB/optimiation/'
matfile_target = os.path.join(folder, 'test.mat')  # Where to save the .mat file to
output_target = os.path.join(folder, 'test.txt')    # Where to save the output log

result = sdpt3glue.sdpt3_solve_problem(problem, sdpt3glue.MATLAB, matfile_target,
                                       output_target=output_target)