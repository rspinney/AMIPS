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

#include "integrator.hpp" 

/*****************************************************************/
/*
    MAIN PARTICLE UPDATE FUNCTION
*/
/*****************************************************************/

/* 
* Updates position and orientation of a given particle
*/
void integrator::updateParticle(const integrationInfo &a_info) noexcept 
{
    updateDirection(a_info);
    updatePosition(a_info);
}

/*****************************************************************/
/*
    DIRECTION UPDATE - Poisson update
*/
/*****************************************************************/

/* 
* Updates orientation and tumble time of a given particle
* - sends new value to new buffer at every time step
* - if the value has changed, generates the next flip time
*/
void integrator::updateDirection(const integrationInfo &a_info) noexcept 
{
    const uint32_t &particle = a_info.particle;
    const bool update = (m_state.simulationClockTime() >= m_state.flipTime(particle));
    const uint8_t currentIdx = m_state.dirIndex(particle);
    const uint8_t nextIdx = static_cast<uint8_t>(update * (1u - currentIdx) + (!update) * currentIdx); // new direction index (branchless)
    m_state.updateVelocityIdx(particle, nextIdx);

    if (update)
    {
        const uint32_t &threadId = a_info.threadId;
        const double timeScale = (nextIdx > 0) ? constants::invGammaMinus : constants::invGammaPlus;
        const double uniformRandNum = m_threadObjects[threadId].rng().uniform();
        const double nextFlipTime  = (m_state.simulationClockTime() - (timeScale) * std::log(1.0 - uniformRandNum));
        m_state.updateFlipTime(particle, nextFlipTime);
    }
}

/*****************************************************************/
/*
    POSITION UPDATE
*/
/*****************************************************************/

/* 
* Updates continuous position of a single particle
* - sends position to the new buffer
* - appends change of bin information to updateList
*/
void integrator::updatePosition(const integrationInfo &a_info) noexcept 
{   
    const uint32_t &particle = a_info.particle;
    const uint32_t &threadId = a_info.threadId;
    const double dW = m_threadObjects[threadId].rng().gaussian() * constants::sqrtDt;  // Gaussian noise
    const double dx =  (constants::propulsionStrength * m_state.velocity(particle) * (1.0 - constants::quorumStrength * force(a_info))) * constants::dt
                + constants::noiseStrength * dW;
    //m_state.addTotalDx(particle,dx);
    const double newPos = PBC(m_state.position(particle) + dx);   // new position with PBC
    m_state.updatePosition(particle, newPos);                //  record new position
    m_state.registerBinUpdate(a_info, m_threadObjects[threadId].updateList(), newPos); // register new position with bin structure
}

/*****************************************************************/
/*
    INTERACTION FUNCTIONS
*/
/*****************************************************************/

/* 
* Returns interacion component of the drift for a given particle
* - result of interaction kernel over all particles
*/
double integrator::force(const integrationInfo &a_info) const noexcept 
{
    const uint32_t &particle = a_info.particle;
    const uint32_t &bin = a_info.bin;

    constexpr double fMult = (1.0 / (2.0 * constants::interactDist)); // contribution per particle
    double force = -fMult; //account for self interaction
    const auto &binList = m_state.getBinNeighbours(bin); //get list of neighbouring bins

    for (const auto adjacentBin : binList) //for each neighbouring bin (including self)
    {
        const auto &set = m_state.getBinParticles(adjacentBin); //get list of particles within neighbouring bin
        for (const auto neighbourParticle : set)
        {
            const double fabsDx_ij = wrapFabsDx(m_state.position(particle) - m_state.position(neighbourParticle));
            force += (fabsDx_ij < constants::interactDist) * fMult;
        }
    }
    
    return force;
}
