/*
Copyright (c) 2012 Jakob Progsch, Václav Zeman

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

Extensions and c++17 compliance: Copyright (C) 2022-2023, Richard Spinney.
 
*/

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>               // std::vector
#include <queue>                // std::queue
#include <memory>               // std::shared_ptr
#include <thread>               // std::thread
#include <mutex>                // std::mutex
#include <condition_variable>   // std::condition_variable
#include <future>               // std::future
#include <functional>           // std::function
#include <stdexcept>            // throw
#include <utility>              // std::move, std::forward

class threadPool 
{
public:
    threadPool(std::size_t);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::invoke_result<F,Args...>::type>;
    // auto enqueue(F&& f, Args&&... args) 
    //     -> std::future<typename std::result_of<F(Args...)>::type>;
    ~threadPool();
    void restart();
    void restart(std::size_t);
private:
    void forceStop();
    void run();
    // need to keep track of threads so we can join them
    std::vector< std::thread > m_workers;
    // the task queue
    std::queue< std::function<void()> > m_tasks;
    // synchronization
    std::mutex m_queueMutex;
    std::condition_variable m_condition;
    size_t m_nThreads;
    bool m_stop,m_kill;
    //new
    
    std::atomic<int> m_nJobsPending;
    std::mutex m_mainMutex;
    std::condition_variable m_mainCondition;
    public:
    void waitUntilCompleted();
};
 
// the constructor just launches some amount of workers
inline threadPool::threadPool(std::size_t a_nThreads) 
    : m_nThreads(a_nThreads), m_stop(false), m_kill(false), m_nJobsPending(0u)
{
    run();
}

inline void threadPool::waitUntilCompleted() 
{
    std::unique_lock<std::mutex> lock(m_mainMutex);
    m_mainCondition.wait(lock);
}

inline void threadPool::run()
{
    for(std::size_t i = 0u; i < m_nThreads; ++i)
        m_workers.emplace_back(
            [this]
            {
                for(;;)
                {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->m_queueMutex);
                        this->m_condition.wait(lock,
                            [this]{ return this->m_kill || this->m_stop || !this->m_tasks.empty(); });
                        if ((this->m_kill)||(this->m_stop && this->m_tasks.empty()))
                            return;
                        task = std::move(this->m_tasks.front());
                        this->m_tasks.pop();
                    }
                    task();
                    if ( --m_nJobsPending <= 0 ) {
                        m_mainCondition.notify_one();
                    }
                }
            }
        );
}

// add new work item to the pool
template<class F, class... Args>
auto threadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::invoke_result<F,Args...>::type>
    //-> std::future<typename std::result_of<F(Args...)>::type>
{
    //new
    m_nJobsPending++;
    //old
    //using return_type = typename std::result_of<F(Args...)>::type;
    using return_type = typename std::invoke_result<F,Args...>::type;
    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        //don't allow enqueueing after stopping the pool
        if(m_stop)
            throw std::runtime_error("enqueue on stopped threadPool");
        else if(m_kill)
            throw std::runtime_error("adding jobs to killed pool.");
        m_tasks.emplace([task](){ (*task)(); });
    }
    m_condition.notify_one();
    return res;
}

// the destructor joins all threads
inline threadPool::~threadPool()
{
    {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        m_stop = true;
    }
    m_condition.notify_all();
    for(std::thread &worker: m_workers)
        worker.join();
}

inline void threadPool::forceStop()
{
    {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        while(!m_tasks.empty())
            m_tasks.pop();
        m_kill = true;
        m_nJobsPending = 0u;
    }
    
    m_condition.notify_all();
    
    for(std::thread &worker: m_workers)
        worker.join();   
}

inline void threadPool::restart()
{
    forceStop();
    m_workers.clear();    
    m_stop = false;
    m_kill = false;
    m_nJobsPending = 0u;
    run();
}

inline void threadPool::restart(std::size_t a_nThreads)
{
    forceStop();
    m_nThreads = std::max(static_cast<std::size_t>(1u), a_nThreads);
    m_workers.clear();    
    m_stop = false;
    m_kill = false;
    m_nJobsPending = 0u;
    run();
}

#endif
