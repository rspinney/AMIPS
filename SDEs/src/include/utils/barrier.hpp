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
    Implements barrier pattern where multiple threads perform (long)
    running functions and synchronise with each other using condition variables
    implemented here with a "wait()" function

    includes task wrapper "generateMultiTaskSynchronisedLoop" to wrap staged synchronisatiom
    around function arguments
*/

#ifndef BARRIER_H
#define BARRIER_H

#include <cassert>              // assert
#include <condition_variable>   // std::condition_variable
#include <mutex>                // std::mutex, std::unique_lock
#include <vector>               // std::vector
#include <array>                // std::array
#include <cstdint>              // uint32_t etc.
#include <iostream>             // cerr etc.

#include "sortpermutation.hpp"

class barrier
{
    // error states
    enum barrierError{noError = 0, stageIdxOOB = 1, tasksNotSet = 2, inconsistentTotalStages = 3};
    
    std::mutex m_mutex;                     // mutex for access to the barrier
    std::condition_variable m_condition;    // cv for blocking/freeing current access thread
    std::array<std::size_t,2> m_tasksLeft;  // (reverse) counting variables - need two to ensure thread safe reset
    std::size_t m_numThreads;               // record of number of threads to synchronise
    int m_countIdx;                         // index for which counting variable we are using
    std::size_t *m_remainingPtr;            // pointer to the counting variable we are using
    barrierError m_error;                   // error flag
    uint32_t m_totalStages;                 // expected number of serial tasks informing expected number of synchronisations
    uint32_t m_generatedTasks;              // counter of number of wrapped tasks generated - needs to match num threads
    
public:

    barrier(const std::size_t a_threads) :
        m_mutex(),
        m_condition(),
        m_tasksLeft({0u, 0u}), // start with an empty (not-running) state - i.e. waiting on 0 tasks
        m_numThreads(a_threads),
        m_countIdx(0u), //start with first counter
        m_error(noError),
        m_totalStages(0u),
        m_generatedTasks(0u) 
    {
        m_remainingPtr = &m_tasksLeft[m_countIdx]; //start with the first counter variable
    }

    // no moving or copying
    barrier(const barrier& barrier) = delete;
    barrier(barrier&& barrier) = delete;
    barrier& operator=(const barrier& a_barrier) = delete;
    barrier& operator=(barrier&& a_barrier) = delete;

    ~barrier() noexcept
    {
        assert(0u == *m_remainingPtr); // on destruction nothing should be waiting
    }

    void reset()
    {
        m_totalStages = 0u;
        m_generatedTasks = 0u;
        m_tasksLeft = {0u,0u};
    }

    void registerTotalStages(const uint32_t a_numStages)
    {
        if (m_generatedTasks > 0) // error if tasks registered before total stages set
            m_error = inconsistentTotalStages;
        m_totalStages = a_numStages;
    }

    bool error() const
    {
        if (m_error)
        {
            printError();
            return true;
        }
        return false;
    }

    void printError() const
    {
        if (m_error == stageIdxOOB)
            std::cerr << "ERROR: could not launch. Tasks with indices greater that total tasks submitted. Exiting." << std::endl;
        else if (m_error == tasksNotSet)
            std::cerr << "ERROR: could not launch. Number of total tasks not registered. Exiting." << std::endl;
        else if (m_error == inconsistentTotalStages)
            std::cerr << "ERROR: could not launch. Total number of tasks changed after registering tasks. Exiting." << std::endl;
    }

    void wait()
	{
		std::unique_lock<std::mutex> lock(m_mutex);
		if (*m_remainingPtr == 0u) //if the counter is at zero it needs resetting
            *m_remainingPtr = m_numThreads; // reset our current count variable.
            
		if (0u == --*m_remainingPtr) //decrement number of tasks remaining, if after decrement =0, then wake up threads
		{
            m_countIdx = 1u - m_countIdx; //swap the index (0 vs 1)
            m_remainingPtr = &m_tasksLeft[m_countIdx]; //swap the counter variable so we dont interfere with any "waits"
			m_condition.notify_all(); //wake up all the threads waiting on the condition variable
		}
		else
		{
            std::size_t *m_localRemainingPtr = m_remainingPtr; //get a pointer to the count variable that is being used (m_remaining will change)
			m_condition.wait(lock, [m_localRemainingPtr]() { return 0u == *m_localRemainingPtr; });// thread stops waiting when tasks remaining hits zero
		}
	}

    // generates a synchronised loop such that each function out of func/idx pair in a_funcIdxPairs (with type void() - i.e. wrap what you need into a lambda)
    // is called with appropriate synchronisation calls such that they occur as staged tasks in the order given by the paired idx out the already registered # of tasks
    // E.g. if we have a_funcIdxPairs = { {f(),1} , {g(),4} }, and m_totalStages = 8 we return a function object which executes

    // while (a_predicate())
    // {   
    //     // in principle some other function is working here...
    //     this->wait();
    //     f();
    //     this->wait();
    //     // in principle some other function is working here...
    //     this->wait();
    //     // in principle some other function is working here...
    //     this->wait();
    //     g();
    //     this->wait();
    //     // in principle some other function is working here...
    //     this->wait();
    //     // in principle some other function is working here...
    //     this->wait();
    //     // in principle some other function is working here...
    //     this->wait();
    // };

    template<typename T, typename U> // T needs to be a vector of void() functions + uint32_t pairs
    auto generateMultiTaskSynchronisedLoop(T &&a_funcIdxPairs, U &&a_predicate)
    {             
        std::vector<uint32_t> taskGaps(a_funcIdxPairs.size());
        const std::vector<std::size_t> perm = sortedPermutation(a_funcIdxPairs,[](auto a, auto b){return a.second < b.second;});
        applyPermutation(a_funcIdxPairs,perm);
        taskGaps.front() = a_funcIdxPairs.front().second;
        for (std::size_t idx = 1u; idx < a_funcIdxPairs.size(); ++idx)
        {
            taskGaps[idx] = a_funcIdxPairs[idx].second - a_funcIdxPairs[idx - 1].second;
            if (a_funcIdxPairs[idx].second >= m_totalStages)
                m_error = stageIdxOOB;
        }
        uint32_t remainingSyncs = m_totalStages - a_funcIdxPairs.back().second;
        if (m_totalStages == 0)
            m_error = tasksNotSet;
        // just capture by value - performance not crucial when building tasks        
        auto wrappedTasks = [=]() noexcept 
        { 
            while (a_predicate()) // continually loop whilst a_predicate returns true
            {   
                for (std::size_t taskIdx = 0u; taskIdx < a_funcIdxPairs.size(); ++taskIdx)
                {
                    for (uint32_t idx = 0u; idx < taskGaps[taskIdx]; ++idx)
                        this->wait(); // wait for previous tasks
                    a_funcIdxPairs[taskIdx].first(); // call the functions
                }
                for (uint32_t idx = 0u; idx < remainingSyncs; ++idx)
                    this->wait();  // wait for subsequent tasks             
            };
        };      
        ++m_generatedTasks;
        m_numThreads = m_generatedTasks;
        return wrappedTasks;
    }
};

#endif /*BARRIER_H*/
