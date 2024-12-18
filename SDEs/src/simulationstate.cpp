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

#include "simulationstate.hpp"

/*
* Initialisation of all data in case of new simulation start
*/
void simulationState::init() noexcept
{
    initialiseParticles();
    initialiseIndexing();
    assignParticlesToThreads();
    m_bins.assignParticlesToBins(m_pos);
}

/*
* Initialisation of all data in case of restarted simulation
*/
void simulationState::setState(const std::vector<double> &a_pos, const std::vector<uint8_t> &a_dirIdxs, const std::vector<double> &a_flipTimes) noexcept
{
    m_pos = a_pos;
    m_dirIdx = a_dirIdxs;
    m_flipTimes = a_flipTimes;
    m_posNext = m_pos; 
    m_dirIdxNext = m_dirIdx;
    initialiseIndexing();
    assignParticlesToThreads();
    m_bins.assignParticlesToBins(m_pos);
}

/*
* Initialises positions, orientations, and tumble times of particles
* - position is distributed either uniformly or in a phase separated state based on parameters
*/
void simulationState::initialiseParticles() noexcept
{
    randNumGen<std::mt19937_64> rng(12345u);

    m_totalDx = std::vector<double>(m_nParticles,0.0L);

    m_pos = std::vector<double>(m_nParticles,0.0L);
    m_dirIdx = std::vector<uint8_t>(m_nParticles);
    m_posNext = std::vector<double>(m_nParticles,0.0L);
    m_dirIdxNext = std::vector<uint8_t>(m_nParticles);
    m_flipTimes = std::vector<double>(m_nParticles);

    for (uint32_t idx = 0u; idx < m_nParticles; ++idx) //generate particle positions
    {
        if (constants::nucleate) // if initial condition is a phase separated state
        {                
            const double highRho = constants::rhoMipsHigh;
            const double lowRho = constants::rhoMipsLow;
            const double highL = (constants::rho - lowRho) / (highRho - lowRho);

            const uint64_t numHigh = static_cast<uint64_t>(((highRho * highL) / (highRho * highL + lowRho * (1.0 - highL))) * static_cast<double>(constants::N));
            const double dxHigh = (constants::L * highL) / static_cast<double>(numHigh);
            const uint64_t numLow = m_nParticles - numHigh ;
            const double dxLow = (constants::L * (1.0 - highL)) / static_cast<double>(numLow);

            if (idx < numHigh)
                m_pos[idx] = 0.5 * (1.0 - highL) * constants::L + idx * dxHigh;
            else if (idx < (numHigh + numLow / 2))
                m_pos[idx] = static_cast<double> (idx - numHigh) *  dxLow;
            else
                m_pos[idx] = static_cast<double> (0.5 * (1.0 + highL) * constants::L + static_cast<double>(idx - numHigh - numLow / 2) *  dxLow);
        }
        else // else:  we homogeneously distribute the particles
        {
            m_pos[idx] = idx * (constants::L / static_cast<double>(m_nParticles));
        }
                    
    }
    for (auto &x : m_dirIdx) // distribute directions according to 2-level system steady state g+ / (g+ + g-)
        x = (rng.uniform() < (constants::gammaPlus  / (constants::gammaPlus + constants::gammaMinus))) ? 1u : 0u;
    
    for (std::size_t idx = 0u; idx < m_flipTimes.size(); ++idx) // distribute times exponentially according to which tumble state they are in
    {
        const double invRate = (m_dirIdx[idx]) ? (1.0 / constants::gammaMinus) : (1.0 / constants::gammaPlus);
        m_flipTimes[idx] = - (invRate) * std::log(1.0 - rng.uniform());
    }

    // initialise swap buffer, mainly for debugging...
    m_posNext = m_pos; 
    m_dirIdxNext = m_dirIdx;

}

/*
* Initialises the record of particle location in buffer
* - at start it is just the identity, i.e. position[i] = i
*/
void simulationState::initialiseIndexing() noexcept
{
    m_particleDataLocations = std::vector<uint32_t>(m_nParticles);
    m_particleIdAtDataLocation = std::vector<uint32_t>(m_nParticles);
    for (uint32_t idx = 0; idx < m_nParticles; ++idx)
    {
        m_particleDataLocations[idx] = idx;
        m_particleIdAtDataLocation[idx] = idx;
    }
}

/*
* Computes instantaneous mean orientation
*/
void simulationState::assignParticlesToThreads() noexcept //initialisation only
{
    m_setBounds = std::vector<std::array<uint32_t,2u> > (m_numThreads);
    const uint32_t binsPerThread = m_bins.numBins() / m_numThreads;
    for (uint32_t idx = 0u; idx < m_numThreads; ++idx)
    {
        m_setBounds[idx][0u] = idx * binsPerThread;
        m_setBounds[idx][1u] = (idx + 1u) * binsPerThread;
    }
    m_setBounds.back()[1u] = m_bins.numBins();
}

/*
* Swaps current data buffer with now completely computed new buffer for next time step
*/
void simulationState::swapBuffers() noexcept
{
    std::swap(m_pos, m_posNext);
    std::swap(m_dirIdx, m_dirIdxNext);
}

/*
* Assigns groups of contiguous bins to threads for computation
* - attempts to keep number of particles per thread approx. constant
*/
void simulationState::updateBinToThreadAssignment() noexcept // assign bins to threads contiguously with roughly equal particles
{
    const uint32_t particlesPerThread = m_nParticles / m_numThreads;
    auto &bounds = m_setBounds;
    const uint32_t numBins = m_bins.numBins();
    uint32_t sumBins = 0u;
    uint32_t threadIdx = 0u;
    uint32_t binIdx = 0u;

    for (uint32_t idx = 0u; idx < numBins; ++idx)
    {
        sumBins += static_cast<uint32_t>(m_bins.getParticles(idx).size());
        if (sumBins >= particlesPerThread)
        {
            sumBins = 0u;
            bounds[threadIdx][0u] = binIdx;
            bounds[threadIdx][1u] = idx;
            binIdx = idx;
            ++threadIdx;
        }
    }
    if (m_numThreads > 1u)
        bounds.back()[0u] = bounds[bounds.size() - 2u][1u];
    else
        bounds.back()[0u] = 0u;
    bounds.back()[1u] =  numBins;
}

/*
* Updates record of particle position in buffer 
* - used following data sorts which change particle ordering
* - only required if outputting trajectories
*/
void simulationState::updateIndexing(const std::vector<std::size_t>& a_permutation) noexcept
{
    std::vector<uint32_t> particleIds = m_particleIdAtDataLocation;
    for (std::size_t idx = 0; idx < a_permutation.size(); ++idx)
    {
        particleIds[idx] = m_particleIdAtDataLocation[a_permutation[idx]];
        m_particleDataLocations[particleIds[idx]] = static_cast<uint32_t>(idx);
    }    
    std::swap(particleIds, m_particleIdAtDataLocation);
}

/*
* Sorts particles according to position and applies the new permutation to all objects
* - improves data locality for performance only
* - reassigns particles to bins based on new index
* - reassigns bins to threads
*/
void simulationState::sortData() // sort the data so that all data for contiguous set of bins is contiguous (i.e. try to achieve data locality)
{    
    const std::vector<std::size_t> sortedIdxs = sortedPermutation(m_pos, std::less<double>());
    if (constants::outputTrajectories)
        m_parallel.submitJob( [&](){this->updateIndexing(std::cref(sortedIdxs));} );
    
    m_parallel.submitJob(&applyPermutation<double>, std::ref(m_pos), std::cref(sortedIdxs));
    m_parallel.submitJob(&applyPermutation<uint8_t>, std::ref(m_dirIdx), std::cref(sortedIdxs));
    m_parallel.submitJob(&applyPermutation<double>, std::ref(m_flipTimes), std::cref(sortedIdxs));
    m_parallel.submitJob(&applyPermutation<double>, std::ref(m_totalDx), std::cref(sortedIdxs));
    m_parallel.waitForJobs();
    m_bins.assignParticlesToBins(m_pos);
    updateBinToThreadAssignment();
}

/*
* Computes mean speed using total accumulated distances
*/
double simulationState::averageSpeed() const noexcept
{
    double sum = 0.0L;
    for (uint32_t idx = 0u; idx < m_nParticles; ++idx)
        sum += totalDx(idx);
    sum /= static_cast<double>(m_nParticles);
    sum /= (m_simulationTime);
    return sum;
}

/*
* Computes instantaneous mean orientation
*/
double simulationState::averageOrientation() const noexcept
{
    double sum = 0.0L;
    for (uint32_t idx = 0u; idx < m_nParticles; ++idx)
        sum += velocity(idx);
    sum /= static_cast<double>(m_nParticles);
    return sum;
}
