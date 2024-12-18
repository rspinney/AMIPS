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

#include "binstructure.hpp"

/* 
* Sets the neighbouring bins for each bin
* - neighbours are bins that need to be searched to compute forces
* - e.g. bin 2 has neighbours 1,2,3
*/
void binStructure::setUpBins() noexcept 
{
    for (uint32_t idx = 0u; idx < m_nBins; ++idx)
    {
        m_binNeighbours[idx][0u] = idx - 1u;
        m_binNeighbours[idx][1u] = idx;
        m_binNeighbours[idx][2u] = idx + 1u;
    }

    // periodic boundaries
    m_binNeighbours.front()[0u] = m_nBins - 1u;
    m_binNeighbours.back()[2u] = 0;
}

/* 
* Initialises the bin record for all particles
* - used at initialisation
*/
void binStructure::assignParticlesToBins(const std::vector<double> &a_pos)  noexcept 
{
    for (auto &x : m_binSets)
        x.clear();
    for (uint32_t idx = 0u; idx < m_nParticles; ++idx)
    {
        uint32_t bin = static_cast<uint32_t>(a_pos[idx] / m_distBin);
        m_binSets[bin].insert(idx);
    }
}

/* 
* Changes bin record for particles associated with given thread
* - changes are not automatically threadsafe
* - use scoped_lock to ensure thread safety
*/
void binStructure::updateBinsThread(std::vector<particleAndBinPair> &updateList) noexcept 
{
    for (auto &x : updateList)
    {
        const auto &particle = x.particle;
        const auto &newBin = x.binNew;
        const auto &oldBin = x.binOld;      
        {
            std::scoped_lock lock(m_mutexes[oldBin], m_mutexes[newBin]); // make operations threadsafe - should automatically avoid deadlock
            m_binSets[oldBin].erase(particle);  // remove from old set
            m_binSets[newBin].insert(particle); // insert into new set
        }
    }
    updateList.clear(); // we have processed the list
}

/* 
* Tests for bin changes and records any changes that need to be implemented for a given particle
* - computes new bin and compares to old bin
* - appends information to updateList if there is a change
*/
void binStructure::registerBinUpdate(const integrationInfo &a_info, std::vector<particleAndBinPair> &updateList, const double &a_posNext) noexcept 
{
    const uint32_t &particle = a_info.particle;
    const uint32_t &oldBin = a_info.bin;
    //register particles with need for bin updates
    uint32_t newBin = static_cast<uint32_t>(a_posNext * m_invDistBin);
    if (newBin != oldBin) //register particle as needing updating
        updateList.emplace_back(particleAndBinPair(particle, oldBin, newBin));
}
