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
/*
    class to encapsulate the state of the simulation including buffers for particle
    position and direction (stored as discrete indices and an array of value), flip 
    times, simulation time.

    conversion of state to a field is delegated to "field"
    linear scaling structures are delegated to "binStructure"
*/

#ifndef PARTICLESTATE_H
#define PARTICLESTATE_H

#include <vector>       // std::vector
#include <array>        // std::array
#include <cstdint>      // uint32_t etc.
#include <cmath>        // std::log
#include <mutex>        // std::mutex, std::scoped_lock
#include <functional>   // std::less<T>

#include "parameters.hpp"
#include "staticasserts.hpp"
#include "random.hpp"
#include "binstructure.hpp"
#include "field.hpp"
#include "parallel.hpp"
#include "sortpermutation.hpp"

class simulationState
{   
    // simulation time   
    const double m_dt;
    double m_simulationTime = 0.0L, m_simulationClockTime = 0.0L;
    uint64_t m_simulationStep = 0UL;
    
    // domain size
    const double m_L;

    // number of particles and threads + thread assignment objects
    const uint32_t m_nParticles, m_numThreads;
    parallel m_parallel; // parallelism for periodically sorting data - responsibility of sim. state
    std::vector<std::array<uint32_t,2u> > m_setBounds; // list of particles in each set - for threads

    // set up functions
    void init() noexcept;
    void initialiseParticles() noexcept;
    void initialiseIndexing() noexcept;
    void assignParticlesToThreads() noexcept;

    // particle state
    std::vector<double> m_pos, m_posNext; // particle positions
    std::vector<double> m_totalDx; // total deviation
    std::vector<uint8_t> m_dirIdx, m_dirIdxNext; // particle direction indidices 0: left, 1: right - for branchless code
    const std::array<double, 2u> m_dirVals = {-1.0L, +1.0L}; // direction values
    std::vector<double> m_flipTimes; // simulation time of next flip
    
    // indexing/tracking particles
    void updateIndexing(const std::vector<std::size_t>&) noexcept;
    std::vector<uint32_t> m_particleDataLocations, m_particleIdAtDataLocation; 
    // particle i is at location j = m_particleDataLocations[i]
    // particle j = particleIdAtDataLocation[i] is at location i

    // linear scaling objects
    binStructure m_bins;

    // singular particle -> continuous field object
    mutable field m_densityField, m_dirField;

    // error status
    bool m_error;

    public:
    //ctor
    simulationState(const uint32_t a_nParticles, const double a_dt, const uint32_t a_numThreads, const double a_L, const double a_maxInteractDist) :
        m_dt(a_dt),
        m_L(a_L),
        m_nParticles(a_nParticles),
        m_numThreads(a_numThreads),
        m_parallel(a_numThreads),
        m_setBounds(a_numThreads),
        m_pos(a_nParticles), m_posNext(a_nParticles),
        m_totalDx(a_nParticles),
        m_dirIdx(a_nParticles), m_dirIdxNext(a_nParticles),
        m_flipTimes(a_nParticles),
        m_bins(a_nParticles,m_L,a_maxInteractDist),
        m_densityField(a_numThreads, constants::fieldResolution, m_L, 1.0 / constants::rhoC),
        m_dirField(a_numThreads, constants::fieldResolution, m_L, 1.0 / constants::rhoC),
        m_error(false)
    {
        init();
        if (m_bins.errorStatus())
            m_error = true;
    }

    bool errorStatus() const noexcept 
    {
        return m_error;
    }

    // retrieve simulation time
    const double& simulationTime() const noexcept{ return m_simulationTime;}
    const double& simulationClockTime() const noexcept{ return m_simulationClockTime;}
    const uint64_t& simulationStep() const noexcept{ return m_simulationStep;}
    //increment simulation time (every loop)
    void incrementTime(){ m_simulationTime += m_dt; m_simulationClockTime += m_dt; ++m_simulationStep;}
    //swap buffers (every loop)
    void swapBuffers() noexcept;
    // poll statistics of particles
    double averageSpeed() const noexcept;
    double averageOrientation() const noexcept;
    // assign bins to threads
    void updateBinToThreadAssignment() noexcept;
    // sort data by position
    void sortData();
    // reinitialise simulation state with read-in data
    void setState(const std::vector<double> &, const std::vector<uint8_t> &, const std::vector<double> &) noexcept;

    uint32_t arrayPosOfParticle(const std::size_t a_idx) const noexcept 
    {
        return m_particleDataLocations[a_idx];
    }

    uint32_t particleIdOfArrayPos(const std::size_t a_idx) const noexcept 
    {
        return m_particleIdAtDataLocation[a_idx];
    }

    void updateBinsThread(std::vector<particleAndBinPair> &a_updateList) noexcept 
    {
        m_bins.updateBinsThread(a_updateList);
    }
    void registerBinUpdate(const integrationInfo &a_info, std::vector<particleAndBinPair> &a_updateList, const double &a_posNext) noexcept 
    {
        m_bins.registerBinUpdate(a_info, a_updateList, a_posNext);
    }
    const auto& densityField(const double a_fieldScale) const noexcept 
    {
        return m_densityField.calcField(m_pos, a_fieldScale);
    }
    const auto& polarisationDensityField(const double a_fieldScale) const noexcept 
    {
        return m_dirField.calcField(m_pos, m_dirIdx, m_dirVals, a_fieldScale);
    }

    // get particle/thread organisation details

    //return std::array of bins which neighbour supplied bin ID
    const auto& getBinNeighbours(const uint32_t a_bin) const noexcept 
    {
        return m_bins.getNeighbours(a_bin);
    }
    //return std::unordered_set of particles in supplied bin ID
    const auto& getBinParticles(const uint32_t a_bin) const noexcept 
    {
        return m_bins.getParticles(a_bin);
    }
    // return std::array containing limits (low,high) of bins assigned to supplied thread ID
    const auto& getSetBounds(const uint32_t a_threadIdx) const noexcept 
    {
        return m_setBounds[a_threadIdx];
    }

    // accessors  "getters"
    const double& position(const uint32_t a_particle) const noexcept 
    {
        return m_pos[a_particle];
    }
    const double& velocity(const uint32_t a_particle) const noexcept 
    {
        return m_dirVals[m_dirIdx[a_particle]];
    }
    const uint8_t& dirIndex(const uint32_t a_particle) const noexcept 
    {
        return m_dirIdx[a_particle];
    }
    const double& flipTime(const uint32_t a_particle) const noexcept 
    {
        return m_flipTimes[a_particle];
    }
    const double& totalDx(const uint32_t a_particle) const noexcept 
    {
        return m_totalDx[a_particle];
    }

    // modifiers "setters"
    void updatePosition(const uint32_t a_particle, const double a_xNew) noexcept 
    {
        m_posNext[a_particle] = a_xNew;
    }
    void updateVelocityIdx(const uint32_t a_particle, const uint8_t a_dirIdxNew) noexcept 
    {
        m_dirIdxNext[a_particle] = a_dirIdxNew;
    }
    void updateFlipTime(const uint32_t a_particle, const double a_flipTimeNew) noexcept 
    {
        m_flipTimes[a_particle] = a_flipTimeNew;
    }
    void addTotalDx(const uint32_t a_particle, const double a_dx) noexcept 
    {
        m_totalDx[a_particle] += a_dx;
    }
    void setSimulationClockTime(double a_time) noexcept
    {
        m_simulationClockTime = a_time;
    }
};

#endif /*PARTICLESTATE_H*/