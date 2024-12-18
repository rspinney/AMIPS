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
    class to encapsulate tracking of particles within spatial "bins"
    to facilitate quasi-linear scaling by allowing interaction only 
    with particles within local regions

    records of particles in bins are stored in hash sets (unordered_set)
    records of bins which particles belong to are stored in a vector
    records of bins other bins need to interact with are stored in a vector of arrays

    when particles move they must call "registerBinUpdate" to make a record
    of any changes of bin the particle is in, which adds to vectors of particle/bin 
    pairs which have changed, stored in separate vectors for independent threads of
    execution

    when all particles have been moved/integrated the hash sets and particle record
    vectors must be updated by calling "updateBinThread" for each thread of execution.

    Each such thread can access any class wide record allowing for data-races, which 
    are prevented with mutexes
*/

#ifndef BINSTRUCTURE_H
#define BINSTRUCTURE_H

#include <vector>           // std::vector
#include <array>            // std::array
#include <unordered_set>    // std::unordered_set
#include <mutex>            // std::mutex
#include <iostream>         // std::cerr
#include <cstdint>          // uint32_t etc.
#include <utility>          // std::move

#include "structures.hpp"

class binStructure
{
    // number of particles to be tracked
    const uint32_t m_nParticles; // number of particles to place in bins

    // bin properties
    const uint32_t m_nBins;      // number of spatial bins
    const double m_distBin;      // size of bins
    const double m_invDistBin;   // inverse of bin size

    // tracking structure objects
    std::vector<std::array<uint32_t,3> > m_binNeighbours; // which bins need to be searched to compute forces
    std::vector<std::unordered_set<uint32_t> > m_binSets; // set of particles in each bin

    // objects for threadsafe access to bins
    std::vector<std::mutex> m_mutexes;

    // error status
    bool m_error;

    void setUpBins() noexcept;

    public:

    binStructure(const uint32_t a_numParticles, const double a_L, const double a_maxInteractDist) :
       m_nParticles(a_numParticles),
       m_nBins(static_cast<uint32_t>(a_L / a_maxInteractDist)),
       m_distBin(a_L / static_cast<double> (m_nBins)),
       m_invDistBin(1.0 / m_distBin),
       m_binNeighbours(m_nBins),
       m_binSets(m_nBins),
       m_mutexes(m_nBins),
       m_error(false)
    {
        if (m_nBins < 4u)
        {
            std::cerr << "Number of bins too small. Interaction distance too large/field too small. Exiting..." << std::endl;
            m_error = true;
            return;
        }
        setUpBins();   
    }
    
    void assignParticlesToBins(const std::vector<double> &) noexcept;
    void updateBinsThread(std::vector<particleAndBinPair> &) noexcept;
    void registerBinUpdate(const integrationInfo &, std::vector<particleAndBinPair> &, const double &) noexcept;
    
    bool errorStatus() const noexcept 
    {
        return m_error;
    }

    //get bin info
    uint32_t numBins() const noexcept 
    {
        return m_nBins;
    }
    
    const auto& getNeighbours(const uint32_t a_bin) const noexcept 
    {
        return m_binNeighbours[a_bin];
    }

    const auto& getParticles (const uint32_t a_bin) const noexcept 
    {
        return m_binSets[a_bin];
    }
};

#endif /*BINSTRUCTURE_H*/
