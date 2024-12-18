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

#include "model.hpp" 

/*****************************************************************/
/*
    MAIN/TOP-LEVEL RUN FUNCTION
*/
/*****************************************************************/

/*
* Main entry point
* - runs in barrier or thread mode depending on parameters
*/
void model::run()
{
    if (m_error)
    {
        std::cerr << "Error. Exiting..." << std::endl;
        return;
    }
    if (constants::barrier) runBarrier();   
    else runThread();
}

/*****************************************************************/
/*
    MAIN INTEGRATION LOOP
    1. update particles
    2. re-organise bin structure
    3. swap buffers
    4. sort data and assign bins to threads to evenly distribute work to threads
        and maintain data locality (every SORTINTERVAL_DEF simulation time)
    5. output data (every OUTPUTSTATUS_DEF and OUTPUTDATA_DEF simulation time)
    6. inrement timestep

    1 & 2 explicitly broken up into tasks handed to individual threads
    3-6 encapsulated as a single threaded "globalTasks()" 
    (which manages its own multithreading albeit with the same 
    parallel object model.m_parallel)

*/
/*****************************************************************/

/*****************************************************************/
/*
    LAUNCH METHODS 
    - regular threading with futures for each task which is a single
      part of an update in a single time step
        or
    - barrier pattern with condition variables such that each task/thread
      continuously loops over all time steps
*/
/*****************************************************************/

/*****************************************************************/
/*
    BARRIER STRUCTURE LAUNCHER
*/
/*****************************************************************/

// uses barrier pattern to synchronise ongoing individual loops, each in its own thread
// idea is to lower overhead of submitting jobs and keep functions acting on given thread tasks on the same core
// we register tasks(s) to be run in parallel along with their synchronisation order out of numTasks global steps, over m_numThreads
// we register tasks(s) to be run in a single thread along with their synchronisation order
// we launch the tasks which internally have been packaged into loops with appropriate synchronisation

/*
* Barrier method entry point
* - registers three main steps
* - launches the threads
*/
void model::runBarrier()
{
    m_parallel.registerTotalStages(3u); //three stages to ongoing loops
    auto predicate = [this]() noexcept {return m_integrator.running();}; // loop condition
    // 1) main parallel functions - integrate/update bins - try to keep on-core data coherent
    // vector of pairs {void(uint32_t i) function, uint32_t} to be run in a single loop, in multiple instances in many threads, with i the thread index
    std::vector<std::pair<std::function<void(const uint32_t)>, uint32_t> > funcIdxsPar; //pair = {function, order index}
    funcIdxsPar.emplace_back(std::make_pair([this](const uint32_t a_threadIdx) noexcept {m_integrator.updateParticleSet(a_threadIdx);}, 0u));
    funcIdxsPar.emplace_back(std::make_pair([this](const uint32_t a_threadIdx) noexcept {m_integrator.updateBinSet(a_threadIdx);}, 1u));
    m_parallel.registerMultiThreadMultiStageLoops(std::move(funcIdxsPar), predicate, m_numThreads);
    // 2) single thread function for remaining tasks
    // vector (here just one function) of pairs {void(uint32_t i) function, uint32_t} to be run in a loop in a single thread
    std::vector<std::pair<std::function<void()>, uint32_t> > funcIdxsSingle = {std::make_pair([this]() noexcept {m_integrator.globalTasks();}, 2u)};//pair = {function, order index}
    m_parallel.registerSingleThreadMultiStageLoop(std::move(funcIdxsSingle), predicate);
    m_parallel.launchRegisteredTasks(); //run the above loops
}

/*****************************************************************/
/*
    REGULAR THREADED ENTRY POINT AND MAIN LOOP
*/
/*****************************************************************/

// simpler interface for reference - tasks separately sent to thread pool on every timestep
// overhead of continually submiting to threadpool (arguably same as barrier)
// but also doesnt't keep all thread data (rngs, etc) in one function executing one loop
// thus expect poorer cache performance etc. etc.
// explicitly - core 1 might get task 1 for integrating, but the same core might get task 2 for
// updating the bins. And core 1 might have the necessary data already in cache. Then same argument
// for reverting back to integration.

/*
* Regular thread method entry point
* - single predicate loop
* - calls main tasks which each manage parallelisation
*/
void model::runThread()
{
    while (m_integrator.running())
    {
        m_integrator.updateParticles(m_parallel); // integrate
        m_integrator.updateBins(m_parallel);      // update tracking structure
        m_integrator.globalTasks();               // swap buffers, output, increment time
    }
}
