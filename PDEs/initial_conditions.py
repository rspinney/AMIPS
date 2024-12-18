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

import numpy as np
import math
from mpi4py import MPI

########################################
########### INITIAL CONDITIONS #########
########################################

# Class providing initial condition data for the system in a piece wise constant phase separated state
class initial_separated_state():

    def __init__(self,rho0,L,beta,rho_liquid,rho_gas,repeats, Pe = 0, v = 0, J = 0):
        self.beta = beta
        self.L = L
        self.rho_gas = rho_gas
        self.rho_liquid = rho_liquid
        self.vol_liquid = (self.rho_gas - rho0) / (self.rho_gas - self.rho_liquid)
        self.repeats = np.maximum(1, np.floor(repeats))
        self.region_L = self.L / self.repeats
        self.Pe = Pe
        self.v = v
        self.J = J

        if MPI.COMM_WORLD.rank == 0:
            print("Two densities: (" + str(self.rho_gas) + ", " + str(self.rho_liquid) + ") requires volume of " + str(self.vol_liquid) + " to achieve average density " + str(rho0))
        if self.vol_liquid <= 0 or self.vol_liquid >= 1:
            raise Exception("Unphysical initial condition requested.") 

    def description(self):
        return "2-Phase piecewise continuous profile formed of " + str(self.repeats) + " repeats of pairs of contiguous regions wth densities " + str(self.rho_gas) + " and " + str(self.rho_liquid)

    # density profile
    def rho_init(self,x):
        values = np.zeros((1, x.shape[1]))
        for i, val in enumerate(x[0]):
            dx = math.fmod(x[0][i] + 0.75 * self.region_L, self.region_L) / self.region_L
            if (dx < self.vol_liquid): 
                values[0][i] = self.rho_liquid
            else:
                values[0][i] = self.rho_gas
        return values
    
    # polarisation density profile for full simulation
    def psi_init(self,x):
        values = np.zeros((1, x.shape[1]))
        for i, val in enumerate(x[0]):
            dx = math.fmod(x[0][i] + 0.75 * self.region_L, self.region_L) / self.region_L
            if (dx < self.vol_liquid):
                values[0][i] = self.rho_liquid * (1 - self.beta) / (1 + self.beta)
            else:
                values[0][i] = self.rho_gas * (1 - self.beta) / (1 + self.beta)
        return values
    
    # mu field profile for implicit NSAMB simulation
    def mu_init(self,x):
        values = np.ones((1, x.shape[1]))
        rho = self.rho_init(x)
        rho = rho[0]
        for i, val in enumerate(x[0]):
            values[0][i] = self.Pe * rho[i] * (1 - rho[i]) - self.v * (self.v * rho[i] + self.J) / (self.Pe * (1 - rho[i])) - ((1 + self.beta) * (1 + self.beta) / (2 * self.beta * self.Pe)) * np.log(1 - rho[i])
        return values

# Class providing initial condition data for the system in constant state with a Gaussian perturbation
class initial_wrapped_gaussian_state():

    def __init__(self,rho0,amplitude,sigma,L,beta,repeats, Pe = 4.0, v = -1.0, J = 1.0):
        self.beta = beta
        self.L = L
        self.sigma = sigma
        self.sigma_input = sigma
        self.amplitude = amplitude #* (self.sigma * np.sqrt(2 * np.pi))
        self.rho0 = rho0
        self.repeats = repeats
        self.Pe = Pe
        self.v = v
        self.J = J
        if self.rho0 + self.rho0 * self.amplitude / (self.sigma * np.sqrt(2 * np.pi)) > 1 or self.rho0 + self.rho0 * self.amplitude / (self.sigma * np.sqrt(2 * np.pi)) < 0:
            raise Exception("Unphysical initial condition requested.") 

    def description(self):
        return "Wrapped gaussian profile with standard deviation = " + str(self.sigma_input) + " with formula rho = (1 - amplitude / L) * rho_0 + amplitude * rho0 * N_wrap(0,sigma^2), with amplitude = " + str(self.amplitude)

    # density profile
    def rho_init(self,x):
        values = self.rho0 * (1 - self.repeats * self.amplitude / self.L) * np.ones((1, x.shape[1]))
        for i, v in enumerate(x[0]):
            for spike in range(0,self.repeats):
                x_spike = (spike + 1/2) * self.L / (self.repeats)
                dx = v - x_spike
                values[0][i] += self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))
                for j in range(1,5):
                    dx = v - x_spike + j * self.L
                    values[0][i] += self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))
                    dx = v - x_spike - j * self.L
                    values[0][i] += self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))       
        return values
    
    # polarisation density profile for full simulation
    def psi_init(self,x):
        values = ((1 - self.beta) / (1 + self.beta)) * self.rho0 * (1 - self.repeats * self.amplitude / self.L) * np.ones((1, x.shape[1]))
        for i, v in enumerate(x[0]):            
            for spike in range(0,self.repeats):
                x_spike = (spike + 1/2) * self.L / (self.repeats)
                dx = v - x_spike
                values[0][i] += ((1 - self.beta) / (1 + self.beta)) * self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))
                for j in range(1,5):
                    dx = v - x_spike + j * self.L
                    values[0][i] += ((1 - self.beta) / (1 + self.beta)) * self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))
                    dx = v - x_spike - j * self.L
                    values[0][i] += ((1 - self.beta) / (1 + self.beta)) * self.amplitude * self.rho0 * np.exp( - (dx * dx) / (2 * self.sigma * self.sigma)) / (self.sigma * np.sqrt(2 * np.pi))              
        return values
    
    # mu field profile for implicit NSAMB simulation
    def mu_init(self,x):
        values = np.ones((1, x.shape[1]))
        rho = self.rho_init(x)
        rho = rho[0]
        for i,val in enumerate(x[0]):            
            values[0][i] = self.Pe * rho[i] * (1 - rho[i]) -self.v*(self.v * rho[i] + self.J) / (self.Pe * (1 - rho[i])) - ((1 + self.beta) * (1 + self.beta) / (2 * self.beta * self.Pe)) * np.log(1 - rho[i])
        return values

