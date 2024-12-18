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
    Wraps barrier and threadpool objects to encapsulate calling and synchronisation
    of mult-threaded tasks
*/

#ifndef PARALLEL_H
#define PARALLEL_H

#include <cstdint>     // uint32_t etc.
#include <queue>       // std::queue
#include <vector>      // std::vector
#include <functional>  // std::function
#include <future>      // std::future
#include <utility>     // std::forward, std::move

#include "threadpool.hpp"
#include "barrier.hpp"

class parallel
{
    uint32_t m_numThreads; //number of threads (working at one time)
    threadPool m_pool;
    barrier m_barrier;
    std::queue<std::function<void()> > m_tasks;
    std::vector<std::future<void> > m_results;

public:
    
    parallel(const uint32_t a_numThreads) : m_numThreads(a_numThreads), m_pool(a_numThreads), m_barrier(a_numThreads){} 
    ~parallel(){}

    uint32_t numThreads() {return m_numThreads;}

    /*********************************************/
    /**************DIRECT INTERFACE***************/
    /*********************************************/
    
    // submit a single function returning a templated future to thread pool
    // user has to manually deal with future
    template<class T, class... Args>
    auto enqueue(T&& f, Args&&... args) 
    -> std::future<typename std::result_of<T(Args...)>::type>
    {
        return m_pool.enqueue(std::forward<T>(f), std::forward<Args>(args)...);
    }

    // wait for other threads calling sychronise
    void synchronise()
    {
        m_barrier.wait();
    }

    /*********************************************/
    /**************INDIVIDUAL TASKS***************/
    /*********************************************/

    // submit individual void jobs straight to the pool - class handles the future
    template<typename T, class... Args>
    void submitJob(T a_func, Args&&... args) 
    {
        m_results.emplace_back(m_pool.enqueue(std::forward<T>(a_func), std::forward<Args>(args)...));
    }

    // wait for all submitted jobs
    void waitForJobs()
    {
        for (auto &res : m_results)
            res.wait();
        m_results.clear();
    }

    /*********************************************/
    /****************INDEXED TASKS****************/
    /*********************************************/

    // run a thread-indexed function in parallel over m_numThreads threads
    // synchronises when all threads are finished
    // i.e. will run a_func(i) from i = 0 to i = numThreads - 1, each in a separate thread
    template<typename T>
    void parallelIndexedTask(T &&a_func)
    {
        std::vector<std::future<void> > results;
        for (uint32_t idx = 0u; idx < m_numThreads; ++idx) // not perfect forwarding here, due to re-use
            results.emplace_back(m_pool.enqueue(a_func, idx));

        for (auto &res : results)
            res.wait();
    }

    /*********************************************/
    /*********LOOP SYNCHRONISATION TASKS**********/
    /*********************************************/

    void reset()
    {
        m_pool.restart(m_numThreads);
        m_barrier.reset();
        while(!m_tasks.empty()) 
            m_tasks.pop();
    }

    void registerTotalStages(const uint32_t a_numStages)
    {
        m_barrier.registerTotalStages(a_numStages);
    }

    // register several thread-indexed functions to run and synchronise in parallel over a_numThreads threads, 
    // to loop until predicate is no longer true
    // synchronising with all other registered functions at steps according to the idxs in second pos in pair in a_funcIdxPairs
    // out of a_totalIndex total such tasks
    template <typename T, typename U> // T is a vector of std::pairs of void(const uint32_t) functions/lambdas and uint32_t
    void registerMultiThreadMultiStageLoops(T && a_funcIdxPairs, U && a_predicate, uint32_t a_numThreads)
    {    
        for (uint32_t idx = 0u; idx < a_numThreads; ++idx)
        {   // convert indexed tasks with type void(int) to un-indexed tasks with type void()
            std::vector<std::pair<std::function<void()>, uint32_t> > tasks;
            for (auto &funcIdxPair : a_funcIdxPairs)
                tasks.emplace_back(std::make_pair([=]() noexcept {funcIdxPair.first(idx);}, funcIdxPair.second));
            auto task = m_barrier.generateMultiTaskSynchronisedLoop(std::move(tasks), a_predicate);
            m_tasks.push(std::move(task));
        }
    }

    // register several functions to run and synchronise in parallel over a single thread, 
    // to loop until predicate is no longer true
    // synchronising with all other registered functions at step equal to the idx in second pos in pair in a_funcIdxPairs
    // out of a preregistered total number of such tasks
    template <typename T, typename U> // T is a vector of std::pairs of void(const uint32_t) functions/lambdas and uint32_t
    void registerSingleThreadMultiStageLoop(T && a_funcIdxPairs, U && a_predicate)
    {    
        auto task = m_barrier.generateMultiTaskSynchronisedLoop(std::forward<T>(a_funcIdxPairs), a_predicate);
        m_tasks.push(std::move(task));
    }

    // launch all registered functions
    void launchRegisteredTasks()
    {        
        const std::size_t numThreadsRequired = m_tasks.size(); // number of waiting threads
        m_pool.restart(numThreadsRequired); //number of total running threads
        if (m_barrier.error())
            return;
        std::vector<std::future<void> > results;
        while (!m_tasks.empty())
        {
            auto task = m_tasks.front();
            results.emplace_back(m_pool.enqueue(std::move(task)));
            m_tasks.pop();
        }
        for (auto &future : results)
            future.wait();
    }
};

#endif /*PARALLEL_H*/
