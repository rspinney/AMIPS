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
    class to encapsulate file handling and output
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include<unistd.h> 

#include <vector>       // std::vector
#include <string>       // std::string
#include <sstream>      // std::stringstream
#include <fstream>      // std::ofstream
#include <iomanip>      // std::setw, std::setprecision
#include <filesystem>   // std::filesystem::create_directories

#include "parameters.hpp"
#include "staticasserts.hpp"
#include "timer.hpp"
#include "simulationstate.hpp"
#include "parallel.hpp"

class io
{
    class file // file class for keeping file state etc.
    {
        std::string m_str; // file path 
        std::fstream m_file; // file stream
        public:
        file(){}
        ~file(){m_file.close();}
        void setPath(const std::string a_str){m_str = a_str;}
        void openForWrite(){m_file.open(m_str.c_str(), std::fstream::out);}
        void openForAppend(){m_file.open(m_str.c_str(), std::fstream::app | std::fstream::out);}
        bool openForRead()
        {
            try 
            {
                m_file.open(m_str.c_str(), std::fstream::in);
            }
            catch(const std::exception& e)
            {
                std::cerr << "Error opening file: " << m_str <<std::endl;
                std::cerr << "Encountered exception: " << e.what() << std::endl;
                return false;
            }
            return true;            
        }
        void close() noexcept {m_file.close();}
        const std::string& path() noexcept {return m_str;}
        std::fstream& stream() noexcept {return m_file;}
        std::fstream& getline(std::string &) noexcept;
    };

    //files and paths
    std::string m_directoryStr, m_fullDirectoryStr, m_timestampStr;
    file m_status, m_trajectoryTime, m_trajectoryPos, m_trajectoryDir, m_parameter, m_densityField, m_dirField, m_leftField, m_rightField, m_state;

    //append flag
    bool m_append;

    //field parameters for output
    const double m_fieldScale = constants::fieldScale; //basis function lengthscale

    //output timing properties
    double m_outputInterval, m_outputStateInterval, m_outputDataInterval, m_outputTrajectoryInterval;
    uint64_t m_outputIntervalSteps, m_outputStateIntervalSteps, m_outputDataIntervalSteps, m_outputTrajectoryIntervalSteps;
    mutable uint64_t m_outputCount = 0u;

    //timer
    mutable timer clock;

    // open files
    void init();
    // query time and write to files
    void outputStatus(const simulationState &) noexcept;
    void outputData(const simulationState &) noexcept;
    void outputTime(const simulationState &, file &) noexcept;
    void outputPositions(const simulationState &, file &) noexcept;
    void outputDirections(const simulationState &, file &) noexcept;
    void outputFlipTimes(const simulationState &, file &) noexcept;
    void outputParticleState(const simulationState &, file &, uint32_t) noexcept;
    void outputSpecies(const std::vector<double> &, const std::vector<double> &) noexcept;
    void outputField(const std::vector<double> &) noexcept;
    void outputPolarisationField(const std::vector<double> &) noexcept;
    // true/false write output funcs
    bool outputStatusCondition(const simulationState &) const noexcept;
    bool outputDataCondition(const simulationState &) const noexcept;
    bool outputStateCondition(const simulationState &) const noexcept;
    bool outputTrajectoriesCondition(const simulationState &) const noexcept;
    // write params to file
    void outputParameters() noexcept;
    // output tasks
    void outputFields(const simulationState &) noexcept;
    void outputTrajectories(const simulationState &) noexcept;
    void outputState(const simulationState &) noexcept;
    
    public:

    void processOutput(const simulationState &a_state) noexcept
    {
        outputStatus(a_state);
        outputData(a_state);
    }

    void initialOutput(const simulationState &a_state) noexcept
    {
        if (constants::burn < constants::dt)
        {
            if (constants::outputState)
                outputState(a_state);
            if (constants::outputField || constants::outputDirField || constants::outputSpeciesField)
                outputFields(a_state);
            if (constants::outputTrajectories)
                outputTrajectories(a_state);
        }
    }

    //read state
    bool readState(simulationState &, const std::string);

    io(const std::string a_directoryStr, const std::string a_storedStateStr) :
        m_directoryStr(a_directoryStr),
        m_outputInterval(constants::outputInterval),
        m_outputStateInterval(constants::outputStateInterval),
        m_outputDataInterval(constants::outputDataInterval),
        m_outputTrajectoryInterval(constants::outputTrajInterval)
    {
        m_outputIntervalSteps = static_cast<uint64_t>(m_outputInterval / constants::dt);
        m_outputStateIntervalSteps = static_cast<uint64_t>(m_outputStateInterval / constants::dt);
        m_outputDataIntervalSteps = static_cast<uint64_t>(m_outputDataInterval / constants::dt);
        m_outputTrajectoryIntervalSteps = static_cast<uint64_t>(m_outputTrajectoryInterval / constants::dt);

        if (a_storedStateStr.size()) // we are reading state and appending
            m_append = true;
        else
            m_append = false;

        init();
        outputParameters();
        clock.start();
    }
    ~io()
    {
        m_status.close();
        m_parameter.close();
        if (constants::outputDirField)
            m_dirField.close();
        if (constants::outputField)
            m_densityField.close();
        if (constants::outputSpeciesField)
        {
            m_leftField.close();
            m_rightField.close();
        }
        if (constants::outputTrajectories)
        {
            m_trajectoryTime.close();
            m_trajectoryPos.close();
            m_trajectoryDir.close();
        }
        clock.stop();
    }
};

#endif /* OUTPUT_H */
