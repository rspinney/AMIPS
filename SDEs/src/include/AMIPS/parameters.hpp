/*****************************************************************************/
/***************** Copyright (C) 2022-2023, Richard Spinney. *****************/
/*****************************************************************************/
//                                                                           //
//    This program is free software: you can redistribute it and/or modify   //
//    it under the terms of the GNU General Public License as published by   //
//    the Free Software Foundation, either version 3 of the License, or      //
//    (at your option) any later version.                                    //
//                                                                           //
//    This program is distributed in the hope that it will be useful,        //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//    GNU General Public License for more details.                           //
//                                                                           //
//    You should have received a copy of the GNU General Public License      //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdint>
#include <cmath>
//#include <string_view>
#include "constfuncs.hpp"

namespace constants
{

/***************************************/
/*************MATH-CONSTANTS************/
/***************************************/
 
constexpr double pi = 3.14159265359;

/***************************************/
/*************MULTI-THREADING***********/
/***************************************/

constexpr uint32_t numThreads = 12; // degree of parallelisation
// synchronisation option:
// barrier object (true) - all threads have one task which periodically waits on condition variables
// threadpool tasks (false) - each thread is only given task up to next synchronisation event
constexpr bool barrier = true;

/***************************************/
/***************DEBUGGING***************/
/***************************************/

// if true, performs an inefficient aggregation and sort of particles in each bin before integration
// to get deterministic ordering of execution (same behaviour each time), necessary 
// due to non-deterministic ordering of std::unordered_set 
// for debugging only
constexpr bool deterministic = false;

/***************************************/
/**************NUM PARTICLES************/
/***************************************/

constexpr uint32_t N = 320000u;

/***************************************/
/***************DOMAIN SIZE*************/
/***************************************/

constexpr double L = 128.0 * constants::pi;

/***************************************/
/**********INTEGRATION PARAMETERS*******/
/***************************************/

constexpr double dt = 0.001;                         /*time step in simulation seconds*/
constexpr double sqrtDt = sqrtConst(constants::dt);  /*square root of time step using compile time sqrt*/

/***************************************/
/*************OUTPUT-OPTIONS************/
/***************************************/

constexpr bool timeStampPath = false;       /* appends data path with current date/time.*/
constexpr bool outputField = true;          /* output density field - output file: field.dat */
constexpr bool outputDirField = false;      /* output polarisation-density field - output file dirfield.dat */
constexpr bool outputSpeciesField = false;  /* output species (left/right) density fields - output files: leftfield.dat, rightfield.dat */
constexpr bool outputTrajectories = false;  /* output individual trajectories - ouput files: trajectoriesTime.dat, trajectoriesPos.dat, trajectoriesDir.dat*/
constexpr bool outputState = true;          /* output (and overwrite) full particle state (for saving simulation state) - output file: recentState.dat */
// output of trajectories requires identification of each particle ID - uses additional bookkeeping at slight penalty
// note  - output of trajectories leads to large files, use caution. Intended for debugging etc.

/***************************************/
/*******TIME-PARAMETERS/SETTINGS********/
/***************************************/

/*****PERFORMANCE******/

constexpr double sortInterval = 0.2; /*time between sorts for data locality*/
constexpr uint64_t sortIntervalSteps = constants::sortInterval / constants::dt; /*time steps between sorts for data locality*/

/****GENERAL OUTPUT****/

constexpr double outputInterval = 1.0;                                  /*time in simulation seconds between status updates*/
constexpr double outputTrajInterval = constants::dt;                    /*time in simulation seconds between trajectory output (to trajectoriesPos.dat etc.)*/
constexpr double outputDataInterval = 1.0;                              /*time in simulation seconds between data output (to field.dat etc.)*/
constexpr double outputStateInterval = constants::outputDataInterval;   /*time in simulation seconds between state output (to recentState.dat)*/

constexpr double burn = 0.0;                    /*time in simulation seconds before data output*/
constexpr double totalTime = 4000000.0;         /*time of whole simulation in seconds INCLUDING burn time*/ 

/***************************************/
/*******DENSITY/INITIAL-CONDITIONS******/
/***************************************/

constexpr double rhoSpec = 0.65;
constexpr bool nucleate = true;      /*start sim in phase separated state true/false*/
constexpr double rhoMipsHigh = 0.99; //0.99  /*liquid phase starting density*/
constexpr double rhoMipsLow = 0.4;  //0.22 /*gas phase starting density*/
constexpr double rho = (constants::nucleate) ? (0.5 * (constants::rhoMipsHigh + constants::rhoMipsLow)) : constants::rhoSpec;

/***************************************/
/*************MODEL PARAMETERS**********/
/***************************************/

/******UNIVERSAL*******/

constexpr double rhoC = (1.0 / (constants::rho * constants::L));
constexpr double quorumStrength = (1.0 / (constants::rhoC * static_cast<double>(constants::N)));

/*****USER-DEFINED*****/

constexpr double propulsionStrength = (6.0);
constexpr double diffusionConst = (1.0);
constexpr double noiseStrength = sqrtConst(2.0 * constants::diffusionConst); //(1.41421356237); /*noiseStrength := sqrt(2*D)*/
constexpr double beta = 0.5;
constexpr double Pe = 7.0;


/***************************************/
/********TUMBLING RATES - RESULTANT*****/
/***************************************/

constexpr double gammaMinus = (((1.0 + constants::beta) * (constants::propulsionStrength * constants::propulsionStrength)) / (constants::noiseStrength * constants::noiseStrength * constants::Pe * constants::Pe));
constexpr double gammaPlus = (constants::gammaMinus / constants::beta);
constexpr double invGammaMinus = (1.0 / (constants::gammaMinus));
constexpr double invGammaPlus = (1.0 / (constants::gammaPlus));

/***************************************/
/**********INTERACTION PARAMETERS*******/
/***************************************/

constexpr double interactDist = 0.005;                               /* distance between particles when quorum sensing starts*/
constexpr double maxInteractDist = (1.01 * constants::interactDist); /* distance from particles checked for other particles when computing forces*/

/***************************************/
/***********OUTPUT/CG PARAMETERS********/
/***************************************/

constexpr double fieldScale = (15.0 * constants::interactDist);                                /* width of gaussian basis functions for field output*/
constexpr uint32_t fieldPixels = 8000;                                                         /* number of data points in the field output*/
constexpr double fieldResolution = ((1.0 / static_cast<double>(fieldPixels)) * constants::L);  /* resolution of empirical field exported as data output*/

}

#endif /*PARAMETERS_H*/
