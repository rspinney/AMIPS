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

from mpi4py import MPI
from dolfinx.fem import Function, functionspace
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.mesh import (create_interval, locate_entities_boundary, meshtags)
from basix.ufl import (element, mixed_element)
import numpy as np
import math
from ufl import (TestFunctions, dx, inner, split, derivative, ln)
import dolfinx_mpc
from nonlinear_assembly import *
from initial_conditions import *
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import time
from enum import Enum

########################################
############### WHICH PDE ##############
########################################

class Theory(Enum):
    Full = 1           # Full dynamic in rho and psi
    ImplicitNSAMB = 2  # single order parameter theory with implicit parameters v and J (set as constants)

########################################
############# MAIN FUNCTION ############
########################################

def run_aMIPS(NSAMB, L, T, mesh_size, Pe, beta, rho0, dt, dir_str, non_dimensional, output_history, initial_condition, data_interval, output_interval, v = -1.0, J = 1.0, D = 1.0, motility = 6.0):

########################################
############### CONSTANTS ##############
########################################

    # MPI info
    comm = MPI.COMM_WORLD

    # tumbling rates from system parameters
    gm = motility * motility * (1 + beta) / (2 * D * Pe * Pe)
    gp = gm / beta

    # domain
    x_left = 0.0
    x_right = L

    # name strings for output depending on model choice
    if NSAMB == Theory.ImplicitNSAMB:
        aux_str = "mu"
    else:
        aux_str = "Psi"

########################################
############### COUNTERS ###############
########################################

    loop_count = 0
    output_step = int(output_interval / dt) # number of timesteps per output interval
    data_step = int(data_interval / dt) # number of timesteps per output interval

########################################
########## BOUNDARY FUNCTIONS ##########
########################################

    # are we on the boundary
    def periodic_boundary(x):
        return np.isclose(x[0], x_right)

    # periodic relation
    def periodic_relation(x):
        out_x = np.zeros(x.shape)
        out_x[0] = x_right - x[0]
        return out_x

########################################
######## FUNCTION SPACE OBJECTS ########
########################################

    # Create mesh and define function space
    msh = create_interval(comm, points=(x_left, x_right), nx=mesh_size, ghost_mode=dolfinx.cpp.mesh.GhostMode.shared_facet)
    P1 = element("Lagrange", msh.basix_cell(), 1)
    V = functionspace(msh, mixed_element([P1, P1]))
    V_0, _ = V.sub(0).collapse()
    num_points_local = V_0.dofmap.index_map.size_local * V_0.dofmap.index_map_bs
    bcs = [] # no other boundary conditions

    # implement periodic boundary condition
    facets = locate_entities_boundary(msh, msh.topology.dim - 1, periodic_boundary)
    arg_sort = np.argsort(facets)
    mt = meshtags(msh, msh.topology.dim - 1, facets[arg_sort], np.full(len(facets), 2, dtype = np.int32))
    mpc = dolfinx_mpc.MultiPointConstraint(V)
    mpc.create_periodic_constraint_topological(V.sub(0), mt, 2, periodic_relation, bcs, 1.0) # constraint on rho
    mpc.create_periodic_constraint_topological(V.sub(1), mt, 2, periodic_relation, bcs, 1.0) # constraint on aux (psi or mu)
    mpc.finalize()

    # set up solution functions
    u = Function(V)  # current solution
    u_prev = Function(V)  # solution from previous converged step

    # Split mixed functions
    rho, aux = split(u)
    rho_prev, aux_prev = split(u_prev)

    # x values of the mesh
    xvals = sorted(msh.geometry.x[:,0])
    dx_mesh = xvals[1] - xvals[0]

    # all non-autonomous variables must be declared as constant on the mesh to avoid recompilation
    t = dolfinx.fem.Constant(msh, 0.0)
    del_t = dolfinx.fem.Constant(msh, dt)
    
########################################
############## OPEN FILES ##############
########################################

    if comm.rank == 0:  # thread index 0 only
        if output_history == True: # if we are outputing field and output data at each data interval
            rho_file = dir_str + "rho.dat"
            aux_file = dir_str + aux_str + ".dat"
            out_file_rho = open(rho_file, "w")
            out_file_aux = open(aux_file, "w")

            out_file = dir_str + "out.dat"
            stdout_file = open(out_file, "w")

        # outputting order parameter information and most recent field data
        order_file = dir_str + "order.dat"
        out_file_order = open(order_file, "w")
        recent_rho_file = dir_str + "recentRho.dat"
        recent_aux_file = dir_str + "recent" + aux_str + ".dat"

    # initial conditions
    init_rho_file = dir_str + "initialRho.dat"
    init_aux_file = dir_str + "initial" + aux_str + ".dat"

########################################
########## INITIAL CONDITIONS ##########
########################################

    # initial condition for rho
    u.sub(0).interpolate(initial_condition.rho_init)
    
    # initial condition for auxillary field
    if NSAMB == Theory.ImplicitNSAMB:
        u.sub(1).interpolate(initial_condition.mu_init)
    else:
        u.sub(1).interpolate(initial_condition.psi_init)

    u.x.scatter_forward()    

########################################
########## OUTPUT FUNCTIONS ############
########################################

    # get field functions from the mesh form
    def get_profiles(): # has to be called from all threads
        # gather data from all the threads
        rho_dat = msh.comm.gather(u.sub(0).collapse().x.array[:num_points_local])      
        aux_dat = msh.comm.gather(u.sub(1).collapse().x.array[:num_points_local])
        points = msh.comm.gather(msh.geometry.x[:num_points_local,0])        

        xrho, xaux, rho_c, aux_c = [], [], [], []
        if comm.rank == 0: # thread index 0 only
            # concatenate data from threads
            rho_c = np.concatenate(rho_dat)
            aux_c = np.concatenate(aux_dat)
            points_x = np.concatenate(points)
            # pair up x positions with field values
            xrho = list(zip(points_x, rho_c))
            xpsi = list(zip(points_x, aux_c))
            # sort by position
            xrho = sorted(xrho, key = lambda x: x[0])
            xaux = sorted(xpsi, key = lambda x: x[0])
      
        return rho_c, aux_c, xrho, xaux

    # write fields to file
    def write_snapshot(xrho, xaux, rho_file, aux_file):
        outfile = open(rho_file, "w")
        for coord, val in xrho:
            print(str(val), file = outfile)
        outfile.close()

        outfile = open(aux_file, "w")
        for coord, val in xaux:
            print(str(val), file = outfile)
        outfile.close()

    # measure of field deviation
    def accumulated_deviation(rho_dat):
        sum_v = 0.0
        for val in rho_dat:
            sum_v += (val - rho0) * (val - rho0) * dx_mesh
        return sum_v

    # write initial conditions to file
    def write_initial_conditions(initial_rho_file, initial_aux_file):
        rho_dat, aux_dat, xrho, xaux = get_profiles() # has to be done for all threads
        if comm.rank == 0:  # thread index 0 only
            write_snapshot(xrho, xaux, initial_rho_file, initial_aux_file)

    def output():
        # get data profiles
        rho_dat, aux_dat, xrho, xaux = get_profiles() # has to be done for all threads
        if comm.rank == 0: # thread index 0 only
            write_snapshot(xrho, xaux, recent_rho_file, recent_aux_file)
            sum_v = accumulated_deviation(rho_dat)
        
            # append field data to open files
            if output_history == True:
                for coord, val in xrho:
                    print(str(val) + ", ", end = '', file = out_file_rho)
                print(file = out_file_rho, flush = True)
                for coord, val in xaux:
                    print(str(val) + ", ", end = '', file = out_file_aux)
                print(file = out_file_aux, flush = True)
            
            rho_max = np.max(rho_dat)
            rho_min = np.min(rho_dat) 
            diff = rho_max - rho_min

            # find size of largest liquid domain

            # define half-way point between liquid and gaseous domains as decision point
            mid = (rho_max + rho_min) / 2.0

            # find start of first domain
            start = 1
            for idx in range(1,len(rho_dat)-1):
                if rho_dat[idx + 1] >= mid and rho_dat[idx] < mid:
                    start = idx
                    break

            # circ-shift data so it starts at beginning of first domain
            rho_dat = np.roll(rho_dat,-start)

            gasorliq = 1 # domain type flag 0: gas, 1: liquid
            domains = [] # array of length of domains
            start = 1    # index of star of current domain
            for idx in range(1,len(rho_dat)-1):
                if gasorliq == 1:
                    if rho_dat[idx + 1] <= mid: # end of liquid domain
                        finish = idx
                        gasorliq = 0 # we have swapped over into a gas domain
                        domains.append(dx_mesh * (finish - start)) # append length of found domain
                else:
                    if rho_dat[idx + 1] >= mid: # end of gas domain
                        start = idx   # end of gas domain is start of liquid domain
                        gasorliq = 1  # we have swapped over into a liquid domain

            largest_domain = np.max(domains) # find largest of all the domains

            # output to file
            print(str(t.value) + ", " + str(Pe) + ", " + str(diff) + ", " + str(rho_max) + ", " +  str(rho_min) + ", " + str(largest_domain) + ", " +  str(sum_v), file = out_file_order, flush = True)

########################################
######### INITIAL OUTPUT WRITE #########
########################################

    # t=0 output
    output() 
    # write initial conditions to file
    write_initial_conditions(init_rho_file, init_aux_file)

########################################
########### SAVE PARAMETERS ############
########################################

    # save a record of all parameters and write cordinates to file
    if comm.rank == 0: # thread index 0 only
        parameter_file = dir_str + "parameters.dat"
        outfile_parameter = open(parameter_file, "w")

        if NSAMB == Theory.ImplicitNSAMB:
            print("model = NSAMB implicit", file = outfile_parameter)
            print("v = " + str(v), file = outfile_parameter)
            print("J = " + str(J), file = outfile_parameter)
        else:
            print("model = Full (rho + psi)", file = outfile_parameter)

        print("L = " + str(L), file = outfile_parameter)
        print("Pe = " + str(Pe), file = outfile_parameter)
        print("rho0 = " + str(rho0), file = outfile_parameter)
        print("beta = " + str(beta), file = outfile_parameter)
        print("dt = " + str(dt), file = outfile_parameter)
        print("dx = " + str(dx_mesh), file = outfile_parameter)
        print("mesh size = " + str(mesh_size), file = outfile_parameter)
        print("data timestep = " + str(data_interval), file = outfile_parameter)
        print("Initial condition = " + initial_condition.description(), file = outfile_parameter)
        if (non_dimensional == True):
            print("Using non-dimensionalised formulation: partial_t rho = - Pe * partial_x (psi * (1 - rho)) + partial_x^2 rho ...", file = outfile_parameter)
        else:
            print("Using formulation without non-dimensionalisation: partial_t rho = - motility * partial_x (psi * (1 - rho)) + D * partial_x^2 rho ...", file = outfile_parameter)
            print("Additional parameters: D = " + str(D) + " , motility = " + str(motility), file = outfile_parameter)
        outfile_parameter.close()

        # write co-ordinates to file
        xvals_file = dir_str + "xvals.dat"
        out_file_xvals = open(xvals_file, "w")

        for x in xvals:
            print(str(x), file = out_file_xvals)
        out_file_xvals.close()

########################################
########### WEAK FORMULATION ###########
########################################   
    
    theta = 0.5 # crank nicholson
    
    if NSAMB == Theory.ImplicitNSAMB:

        qr, qa = TestFunctions(V)

        mobility_prefactor = (2 * beta * Pe) / ((1 + beta) * (1 + beta))

        rho_mid = (1.0 - theta) * rho_prev + theta * rho

        Jbulk_prev =  Pe * ((1 - beta) / (1 + beta)) * rho_prev * (1 - rho_prev)
        Jbulk = Pe * ((1 - beta) / (1 + beta)) * rho * (1 - rho)
        Jbulk_mid = (1 - theta) * Jbulk_prev + theta * Jbulk

        dmu_times_mobility_mid = mobility_prefactor * ((1 - theta) * (1 - rho_prev) * (aux_prev.dx(0)) + theta * (1 - rho) * (aux.dx(0)))

        # dynamic in rho in terms of bulk current and chemical potential (aux) i.e.
        # \dot{\rho} = -\partial_x ( J_bulk - M \partial_x \mu)
        F0 = -inner (rho - rho_prev, qr) / dt \
        - inner (dmu_times_mobility_mid, (qr).dx(0)) \
        - inner(Jbulk_mid.dx(0),qr) \

        # condition for chemical potential (aux = \mu) i.e.
        # \mu =  Pe*rho*(1-rho)-... i.e. the gradient part of L[\phi] from SI
        F1 = -aux * qa \
            + Pe * rho * (1 - rho) * qa \
            - ((1 + beta) * (1 + beta) / (2 * beta * Pe)) * ln(1 - rho) * qa \
            - v * inner( (rho.dx(0) + v * rho + J) / (Pe * (1 - rho)), qa) \
            + inner( (rho.dx(0) + v * rho + J) / (Pe * (1 - rho)), qa.dx(0))            
            
        F = (F0 + F1) * dx

    else: # Full

        qr, qa = TestFunctions(V)

        if non_dimensional == False:
            advection_constant = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(motility))
            diffusion_constant = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(D))
            reaction_constant1 = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(gm + gp))
            reaction_constant2 = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(gm - gp))
        else:
            advection_constant = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(Pe))
            diffusion_constant = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(1))
            reaction_constant1 = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type( + ((1 + beta) * (1 + beta) / (2 * beta))))
            reaction_constant2 = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type( - ((1 - beta * beta) / (2 * beta))))

        G = dolfinx.fem.Constant(msh, dolfinx.default_scalar_type(0))

        rho_mid = (1.0 - theta) * rho_prev + theta * rho
        aux_mid = (1.0 - theta) * aux_prev + theta * aux
        quorum_aux_mid = (1.0 - theta) * rho_prev * (1 - rho_prev) + theta * rho * (1 - rho)
        quorum_rho_mid = (1.0 - theta) * aux_prev * (1 - rho_prev) + theta * aux * (1 - rho)

        F = \
        ( \
        # time derivative
        inner (rho - rho_prev, qr) / dt \
        + inner (aux - aux_prev, qa) / dt \
        # diffusion terms
        + diffusion_constant * inner ((rho_mid).dx(0), (qr).dx(0)) \
        + diffusion_constant * inner ((aux_mid).dx(0), (qa).dx(0)) \
        # driving terms
        + advection_constant * inner ((quorum_rho_mid).dx(0), qr) \
        + advection_constant * inner ((quorum_aux_mid).dx(0), qa) \
        # tumble terms
        + reaction_constant1 * inner(aux_mid, qa) + reaction_constant2 * inner(rho_mid, qa) \
        #
        # Galilean shift
        - G * inner((rho_mid).dx(0), qr)\
        - G * inner((aux_mid).dx(0), qa)\
        #
        ) * dx

    J = derivative(F, u)
    
    jit_options = {"cffi_extra_compile_args": ["-Ofast", "-march=native"], "cffi_libraries": ["m"]}

    F_compiled = dolfinx.fem.form(F, jit_options=jit_options)
    J_compiled = dolfinx.fem.form(J, jit_options=jit_options)

########################################
########### CONSTRUCT PROBLEM ##########
########################################

    problem = NonlinearMPCProblem(F_compiled, u, mpc, bcs = bcs, J = J_compiled, jit_options = jit_options)
    solver = NewtonSolverMPC(msh.comm, problem, mpc)

    solver.rtol = 1e-7 

    u_prev.x.array[:] = u.x.array
       
########################################
######## PLOT INITIAL CONDITIONS #######
########################################
    
    # # note: doesn't work if running in MPI
    # if comm.rank == 0:
    #     # # examine initial conditions - close to continue
    #     rho_dat = u.sub(0).collapse().x.array
    #     psi_dat = u.sub(1).collapse().x.array
    #     fig1 = plt.figure("Figure 1")
    #     plt.plot(rho_dat)
    #     fig2 = plt.figure("Figure 2")
    #     plt.plot(psi_dat)
    #     # a_dat = u.sub(2).collapse().x.array
    #     # fig3 = plt.figure("Figure 3")
    #     # plt.plot(a_dat)
    #     plt.show()

########################################
########## INTEGRATE/TIME LOOP #########
########################################

    # keep track of time taken
    time_origin = time.time()
    time1 = time_origin

    # time loop
    while (t.value < T):
        # update counters/time
        t.value += del_t.value
        loop_count += 1
        # solve the PDE at this timestep
        r = solver.solve(u)

        # output from thread index 0
        if comm.rank == 0: 
            if (loop_count % output_step == 0):   
                time2 = time.time()        
                out_str = 'on time ' + str(t.value) + '. Speed (Step) = ' + str((output_step * del_t.value) / (time2 - time1)) + ' (unit time / sec)' + '. Speed (Av.) = ' + str(t.value / (time2 - time_origin)) + ' (unit time / sec)'
                print(out_str, flush = True)
                if output_history == True: 
                    print(out_str, file = stdout_file, flush = True)
                time1 = time2
            
        if (loop_count % data_step == 0):
            output()

        # copy data to array for previous values
        u_prev.x.array[:] = u.x.array

    time_final = time.time()

    if comm.rank == 0:
        print("Time taken = " + str(time_final - time_origin) + " seconds.", flush = True)

########################################
################ TIDY UP ###############
########################################

    # close files
    if comm.rank == 0:    
        if output_history == True:
            out_file_rho.close()
            out_file_aux.close()
        out_file_order.close()
