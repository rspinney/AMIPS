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

#ifndef STATICASSERTS_H
#define STATICASSERTS_H

#include <cassert>
#include "parameters.hpp"

static_assert(constants::numThreads > 0u, "Positive number of threads required. Check parameters.\n");
static_assert(constants::N > 0u, "Positive number of particles required. Check parameters.\n");
static_assert(constants::L > 0.0, "Positive domain size required. Check parameters.\n");
static_assert(constants::dt > 0.0, "Positive time step required. Check parameters.\n");
static_assert(constants::interactDist > 0.0, "Positive interaction distance required. Check parameters.\n");
static_assert(constants::maxInteractDist >= constants::interactDist, "Maximum interaction distance must be greater or equal to interaction distance. Check parameters.\n");
static_assert((constants::rho > 0.0) && (constants::rho < 1.0), "Model density must lie in (0,1). Check parameters.\n");
static_assert((constants::rhoMipsHigh >= 0.0) && (constants::rhoMipsHigh <= 1.0), "Initial condition model densities must lie in [0,1]. Check parameters.\n");
static_assert((constants::rhoMipsLow >= 0.0) && (constants::rhoMipsLow <= 1.0), "Initial condition model densities must lie in [0,1]. Check parameters.\n");
static_assert(constants::sortInterval > 0u, "Sorting interval must be positive. Check parameters.\n");
static_assert(constants::outputInterval > 0.0, "Output interval must be positive. Check parameters.\n");
static_assert(constants::outputTrajInterval > 0u, "Trajectory output interval must be positive. Check parameters.\n");
static_assert(constants::outputDataInterval > 0u, "Data output interval must be positive. Check parameters.\n");
static_assert(constants::burn >= 0.0, "Burn interval must not be negative. Check parameters.\n");
static_assert(constants::totalTime > 0.0, "Total simulation time must be positive. Check parameters.\n");
static_assert(constants::propulsionStrength >= 0.0, "Propulsion strength must not be negative. Check parameters.\n");
static_assert(constants::noiseStrength >= 0.0, "Noise strength must not be negative. Check parameters.\n");
static_assert((constants::beta > 0.0) && (constants::beta <= 1.0), "Beta must lie in (0,1]. Check parameters.\n");
static_assert(constants::Pe >= 0.0, "Peclet number must not be negative. Check parameters.\n");
static_assert(constants::gammaPlus >= 0.0, "Left to right tumble rate must not be negative. Check parameters.\n");
static_assert(constants::gammaMinus >= 0.0, "Right to left tumble rate must not be negative. Check parameters.\n");
static_assert(constants::fieldScale > 0.0, "Field scale (basis function width) must be positive. Check parameters.\n");
static_assert(constants::fieldResolution > 0.0, "Field resolution must be positive. Check parameters.\n");

#endif
