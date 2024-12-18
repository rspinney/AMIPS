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
    class to encapsulate the simulation integrator and multithreading
    
    responsibility is to run the main integration loop using methods 
    from integrator using m_parallel
*/

#ifndef MODEL_H
#define MODEL_H

#include <string> // std::string

#include "parameters.hpp"
#include "staticasserts.hpp"
#include "integrator.hpp"
#include "parallel.hpp"

class model
{
    const uint32_t m_numThreads;   // degree of intended parallelism
    integrator m_integrator;       // integration state + methods + output
    parallel m_parallel;           // thread objects/synchronisation
    bool m_error;                  // error flag

    // launchers:
    void runBarrier();  // barrier version
    void runThread();   // thread version
    
public:

    model(const std::string a_folderStr, const std::string a_storedStateStr) :
        m_numThreads(constants::numThreads),
        m_integrator(m_numThreads, constants::N, constants::dt, constants::L, constants::maxInteractDist, a_folderStr, a_storedStateStr),
        m_parallel(m_numThreads),
        m_error(false)
    {
        m_error = m_integrator.errorStatus();
    }
    ~model(){}
    void run(); // entry point
};
#endif /*MODEL_H*/
