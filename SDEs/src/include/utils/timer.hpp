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

#ifndef TIMER_H
#define TIMER_H

/*
    Class wrapper for c++11 timing functionality provided by std::chrono
    Implements simple intervals, total time etc.
*/

#include <chrono>    // std::chrono
#include <mutex>     // std::mutex, std::scoped_lock

class timer 
{
    using timePoint = std::chrono::high_resolution_clock::time_point;
    using tMicroSeconds = std::chrono::microseconds;
    bool m_running, m_stopped; // status flags
    timePoint m_startTime, m_startTimeTemp, m_stopTime, m_stopTimeTemp; // time points
    tMicroSeconds m_mostRecentDurationUsecChrono; // duration in micro seconds
    double m_mostRecentDurationSec; // time in seconds
    std::mutex m_mutex; // mutex for thread safety

public:

    timer() : m_running(false), m_stopped(false){}

    void start() 
    {
        std::scoped_lock lock(m_mutex);
        if (!m_running)
        {
            m_running = true;
            m_startTime = std::chrono::high_resolution_clock::now();
            m_startTimeTemp = m_startTime;
        }
    };

    void stop() 
    {
        std::scoped_lock lock(m_mutex);
        if (m_running) 
        {
            m_stopTime = std::chrono::high_resolution_clock::now();
            m_stopped = true;
        }
    };

    double interval() 
    {
        std::scoped_lock lock(m_mutex);
        m_stopTimeTemp = std::chrono::high_resolution_clock::now();
        if (!m_running) return -1.0;
        m_mostRecentDurationUsecChrono = std::chrono::duration_cast<std::chrono::microseconds>(m_stopTimeTemp - m_startTimeTemp);
        m_mostRecentDurationSec = static_cast<double>(m_mostRecentDurationUsecChrono.count()) / 1000000.0;
        m_startTimeTemp = m_stopTimeTemp;
        return m_mostRecentDurationSec;
    }

    double durationNow() 
    {
        std::scoped_lock lock(m_mutex);
        m_stopTimeTemp = std::chrono::high_resolution_clock::now();
        if (!m_running) return -1.0;
        m_mostRecentDurationUsecChrono = std::chrono::duration_cast<std::chrono::microseconds>(m_stopTimeTemp - m_startTime);
        m_mostRecentDurationSec = static_cast<double>(m_mostRecentDurationUsecChrono.count()) / 1000000.0;
        return m_mostRecentDurationSec;
    }

    double duration() 
    {
        std::scoped_lock lock(m_mutex);
        if (!m_stopped) durationNow();
        m_mostRecentDurationUsecChrono = std::chrono::duration_cast<std::chrono::microseconds>(m_stopTime - m_startTime);
        m_mostRecentDurationSec = static_cast<double>(m_mostRecentDurationUsecChrono.count()) / 1000000.0;
        return m_mostRecentDurationSec;
    };
};

#endif
