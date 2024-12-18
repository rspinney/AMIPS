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
    class to encapsulate integration, 
    
    members m_state and m_io delegate simulation state and output

    parallelism achieved by the interface provided by "parallel.hpp"
*/

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <cstdint>      // uint32_t etc.
#include <vector>       // std::vector
#include <cmath>        // std::fabs, std::log

#include "simulationstate.hpp"
#include "io.hpp"
#include "random.hpp"
#include "structures.hpp"

class integrator
{
    simulationState m_state; // particle state buffers + tracking objects
    io m_io; // input/output structures + methods
    const uint32_t m_numThreads; // degree of parallelism 
    std::vector<threadObjects> m_threadObjects; // thread resources: rngs, particle/bin bounds per thread, bin update lists
    bool m_error; // error flag

    //integration funcs
    void updateParticle(const integrationInfo&) noexcept;
    void updateDirection(const integrationInfo&) noexcept;
    void updatePosition(const integrationInfo&) noexcept;
    double force(const integrationInfo&) const noexcept;

    //inline utility funcs
    double wrapFabsDx(double a_dx) const noexcept
    {
        a_dx = std::fabs(a_dx);
        return (a_dx > 0.5L * constants::L) ? constants::L - a_dx : a_dx;
    }

    double PBC(const double a_x) const noexcept //multiple wrap around version 
    {
        return std::fmod(std::fmod(a_x, constants::L) + constants::L, constants::L);
    } 

    public:

    integrator(const uint32_t a_numThreads, const uint32_t a_nParticles, const double a_dt, const double a_L, const double a_maxDistInteract, const std::string a_folderStr, const std::string a_storedStateStr) :
        m_state(a_nParticles, a_dt, a_numThreads, a_L, a_maxDistInteract),
        m_io(a_folderStr, a_storedStateStr),
        m_numThreads(a_numThreads),
        m_threadObjects(0u),
        m_error(false)
    {
        m_threadObjects.clear();
        for (uint32_t idx = 0; idx < m_numThreads; ++idx)
            m_threadObjects.push_back(threadObjects(idx + 1234u));
        updateSetBounds();
        if (a_storedStateStr.size()) // if we have a non-empty file handle
        {
            if(!m_io.readState(m_state,a_storedStateStr)) // read in old state
                m_error = true;
        }
        else
        {
            m_io.initialOutput(m_state);
        }
        if (m_state.errorStatus())
            m_error = true;
    }

    bool errorStatus() const noexcept
    {
        return m_error;
    }

    void updateSetBounds() noexcept
    {
        for (uint32_t idx = 0; idx < m_numThreads; ++idx)
            m_threadObjects[idx].setBounds() = m_state.getSetBounds(idx);
    }

    /********** MAIN LOOP FUNCS - START **************/

    // main parallel functions - multiple thread callers - STEPS 1 & 2

    // STEP 1: integrate particle
    void updateParticleSet(const uint32_t a_threadIdx) noexcept
    {
        const auto &bounds = m_threadObjects[a_threadIdx].setBounds();
        for (uint32_t bin = bounds[0u]; bin < bounds[1u]; ++bin)
        {
            if (constants::deterministic) // if determinism (debugging) required, slower
            {
                std::vector<uint32_t> particles;
                for (auto particle : m_state.getBinParticles(bin))
                    particles.push_back(particle);
                std::sort(particles.begin(), particles.end());
                for (auto particle : particles)
                    updateParticle({particle, a_threadIdx, bin});
            }
            else // non-determinsitic order, faster
            {
                for (auto particle : m_state.getBinParticles(bin)) // note: order of particles in set is not deterministic
                    updateParticle({particle, a_threadIdx, bin});            
            }          
        }       
    }

    // STEP 2: update bin assignments
    void updateBinSet(const uint32_t a_threadIdx) noexcept
    {
        m_state.updateBinsThread(m_threadObjects[a_threadIdx].updateList());
    }

    // remaining global tasks - single thread caller - STEP 3

    // STEP 3: single thread tasks + rare multithread tasks (e.g. data output, data reshuffling)
    void globalTasks() noexcept
    {
        m_state.swapBuffers();
        if (constants::numThreads > 1u) // no need if the simulation isn't parallel
        {
            if (!(m_state.simulationStep() % constants::sortIntervalSteps))
            {
                m_state.sortData(); //launches its own threads
                updateSetBounds();
            }
        }
        m_io.processOutput(m_state); //launches its own threads
        m_state.incrementTime();
    }

    /********** MAIN LOOP FUNCS - END **************/

    // run main parallel functions over several threads
    void updateParticles(parallel &a_parallel) noexcept
    {
        a_parallel.parallelIndexedTask([this](uint32_t a_threadIdx){updateParticleSet(a_threadIdx);});
    }
    void updateBins(parallel &a_parallel) noexcept
    {
        a_parallel.parallelIndexedTask([this](uint32_t a_threadIdx){updateBinSet(a_threadIdx);});
    }

    //time getter functions
    const double& simulationTime() const noexcept
    {
        return m_state.simulationTime();
    }

    bool running() const noexcept
    {
        return (m_state.simulationTime() < constants::totalTime);
    }
};

#endif /*INTEGRATOR_H*/
