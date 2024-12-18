##############################################################################
################## Copyright (C) 2023-2024, Richard Spinney. #################
##############################################################################
#                                                                            #
#                                                                            #
#    This program is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    This program is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
#                                                                            #
#                                                                            #
##############################################################################

import sys
import os
from mips_fx import *

# print(dolfinx.__version__)

#########################
######## MPI info #######
#########################

comm = MPI.COMM_WORLD

#########################
######### input #########
#########################

num_inputs = len(sys.argv)

# input info/help

if num_inputs > 1:
    if sys.argv[1] == "--help":
        if comm.rank == 0: # if thread 0
            print("\n1D deterministic PDE simulation for non-stationary run and tumble active matter. (c) Richard Spinney 2024.\n\n")
            print("Required input: three floats corresponding to, in order, beta, Pe, rho0. \nFor example `./python3 ./run_amips.py 0.5 7.0 0.7\n")
            print("By default, output written to `./data/$beta_$Pe_$rho0/'\n")
            print("Add suffix to output directory with `-m' flag.\nE.g. `./run_amips.py 0.5 7.0 0.7 -m macroscopic'\nThis will write output to `./data/$beta_$Pe_$rho0_macroscopic/'\n")
        sys.exit()

if num_inputs < 4:
    if comm.rank == 0: # if thread 0
        print("Insufficient arguments. Run `python3 " + os.path.basename(__file__) + " --help' for input instructions\n")
    exit(1)


append_str = ""

# detect optional flag/input for directory suffix

i = 4
while i < num_inputs:
    if sys.argv[i] == "-m":
        if num_inputs < i + 2:
            if comm.rank == 0:
                print("Insufficient arguments. Run `python3 " + os.path.basename(__file__) + " --help' for input instructions\n")
            sys.exit()
        append_str = "_" + sys.argv[i+1]
        i = i + 2
        continue

    if comm.rank == 0:
        print("Unrecognised input. Run `python3 " + os.path.basename(__file__) + " --help' for input instructions\n")
    sys.exit()

append_str += "/"

# strings for the main 3 control parameters
beta_str = sys.argv[1]
Pe_str = sys.argv[2]
rho0_str = sys.argv[3]

#########################
######  parameters ######
#########################

# convert control parameters to numbers
beta = float(beta_str)
Pe = float(Pe_str)
rho0 = float(rho0_str)
 
# model type specification
model = Theory.ImplicitNSAMB # Can take values Theory.Full or Theory.ImplicitNSAMB

# If model == Theory.Full
#   simulates the full system in rho, psi
# If model == Theory.ImplicitNSAMB
#   simulates the single order parameter system dependent on constants v and J which must be specified as parameters

# parameters
L =  4 * 64 * np.pi    # system size
mesh_size = 30000      # number of mesh points
T = 13000              # time horizon
dt = 0.0005            # time step length
data_interval = 1.0    # time between field outputs
output_interval = 1.0  # time between status outputs
output_history = True  # output status updates to out.dat in data directory

# constants required for NSAMB - they are a no-op in the case of the full model
v = -0.83
J = 0.9339

pGsol = 0.495662
pLsol = 0.8820501

v = Pe * ((1 - beta) / (1 + beta)) * (1 - pGsol - pLsol)
J = Pe * ((1 - beta) / (1 + beta)) * (pGsol * pLsol)

#########################
##  initial conditions ##
#########################

#amplitude ref 0.751988
initial_condition = initial_wrapped_gaussian_state(rho0 = rho0, amplitude =  0.001, sigma = 1, L = L, beta = beta, repeats = 1, Pe = Pe, v = v, J = J)
#initial_condition = initial_separated_state(rho0 = rho0, L = L, beta = beta, rhoLiquid = 0.99748, rhoGas = 0.219718, repeats = 1, Pe = Pe, v = v, J = J)

##############################
##  set-up output directory ##
##############################

# construct directory string from control parameters and model type
dir_str = "data/" + "beta_" + beta_str + "_Pe_" + Pe_str + "_p_" + rho0_str
if model == Theory.Full:
    dir_str += "_full"
else:
    dir_str += "_NSAMB"
dir_str += append_str

# write new directory for output
if comm.rank == 0:
    os.makedirs(os.path.dirname(dir_str), exist_ok = True)

#######################
##  main entry point ##
#######################

run_aMIPS( \
NSAMB = model, \
L = L, \
T = T, \
mesh_size = mesh_size,\
Pe = Pe, \
beta = beta, \
rho0 =rho0, \
dt = dt, \
dir_str = dir_str, \
non_dimensional = True, \
initial_condition = initial_condition, \
output_history = output_history, \
data_interval = data_interval, \
output_interval = output_interval, \
v = v, \
J = J, \
D = 1.0 \
)

sys.exit()