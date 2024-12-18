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

#include "io.hpp"

/*
* Reads a line from the file and assigns it to the supplied std::string
* - considers \n, \r, and \r\n all to be end line characters
*/
std::fstream& io::file::getline(std::string& a_str) noexcept 
{
    a_str.clear();    
    //consume string until a new line (of either windows/mac/linux type) 
    std::istream::sentry se(m_file, true);
    std::streambuf* sb = m_file.rdbuf();
    for(;;) 
    {
        int c = sb->sbumpc();
        switch (c) 
        {
            case '\n':
                return m_file;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return m_file;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(a_str.empty())
                    m_file.setstate(std::ios::eofbit);
                return m_file;
            default:
                a_str += static_cast<char>(c);
        }
    }
}

/*
* Creates files for output
* - creates directories if they don't exist.
*/
void io::init()
{
    std::stringstream ss;
    time_t rawtime;
    struct tm *timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);    
    ss << timeinfo->tm_year+1900u << "-" << std::setw(2u) << std::setfill('0') << 
        timeinfo->tm_mon+1u << "-" << std::setw(2u) << std::setfill('0') << timeinfo->tm_mday << 
        " " << std::setw(2u) << std::setfill('0') << timeinfo->tm_hour << 
        ":" << std::setw(2u) << std::setfill('0') << timeinfo->tm_min<< 
        ":" << std::setw(2u) << std::setfill('0') << timeinfo->tm_sec;
    m_timestampStr = ss.str();

    if (constants::timeStampPath)
    {
        m_directoryStr += "_";
        m_directoryStr += m_timestampStr;
        std::replace_if(std::begin(m_directoryStr), std::end(m_directoryStr),
            [](std::string::value_type character) { return character ==' ';},'_');
    }
 
    ss.str("");
    std::string pathStr = std::filesystem::current_path();
    ss << pathStr << "/data/" << m_directoryStr;
    m_fullDirectoryStr = ss.str();

    ss.str("");
    ss << "data/" << m_directoryStr;
    std::filesystem::create_directories(ss.str().c_str());    

    ss.str("");
    ss << "data/" << m_directoryStr << "/parameters.dat";
    m_parameter.setPath(ss.str());

    ss.str("");
    ss << "data/" << m_directoryStr << "/out.dat";
    m_status.setPath(ss.str());
    if (m_append)
        m_status.openForAppend();
    else
        m_status.openForWrite();

    std::cout << "Starting simulation at time: " << m_timestampStr << std::endl;
    std::cout << "Creating data directory: " << m_fullDirectoryStr << std::endl;

    m_status.stream() << "Starting simulation at time: " << m_timestampStr << std::endl;
    m_status.stream() << "Creating data directory: " << m_fullDirectoryStr << std::endl;
    

    if (constants::outputField)
    {
        ss.str("");
        ss << "data/" << m_directoryStr << "/field.dat";
        m_densityField.setPath(ss.str());
        if (m_append)
            m_densityField.openForAppend();
        else
            m_densityField.openForWrite();
        std::cout << "Writing density field to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing density field to file: " << ss.str() << std::endl;
    }

    if (constants::outputDirField)
    {
        ss.str("");
        ss << "data/"<< m_directoryStr << "/dirfield.dat";
        m_dirField.setPath(ss.str());
        if (m_append)
            m_dirField.openForAppend();
        else
            m_dirField.openForWrite();
        std::cout << "Writing polarisation density field to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing polarisation density field to file: " << ss.str() << std::endl;
    }

    if (constants::outputSpeciesField)
    {
        ss.str("");
        ss << "data/" << m_directoryStr << "/leftfield.dat";
        m_leftField.setPath(ss.str());
        if (m_append)
            m_leftField.openForAppend();
        else
            m_leftField.openForWrite();
        std::cout << "Writing left species density field to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing left species density field to file: " << ss.str() << std::endl;

        ss.str("");
        ss << "data/" << m_directoryStr << "/rightfield.dat";
        m_rightField.setPath(ss.str());
        if (m_append)
            m_rightField.openForAppend();
        else
            m_rightField.openForWrite();
        std::cout << "Writing right species density field to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing right species density field to file: " << ss.str() << std::endl;
    }

    if (constants::outputTrajectories)
    {
        ss.str("");
        ss << "data/" << m_directoryStr << "/trajectoriesTime.dat";
        m_trajectoryTime.setPath(ss.str());
        if (m_append)
            m_trajectoryTime.openForAppend();
        else
            m_trajectoryTime.openForWrite();
        std::cout << "Writing trajectory timestamps to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing trajectory timestamps to file: " << ss.str() << std::endl;

        ss.str("");
        ss << "data/" << m_directoryStr << "/trajectoriesPos.dat";
        m_trajectoryPos.setPath(ss.str());
        if (m_append)
            m_trajectoryPos.openForAppend();
        else
            m_trajectoryPos.openForWrite();
        std::cout << "Writing trajectory positions to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing trajectory positions to file: " << ss.str() << std::endl;

        ss.str("");
        ss << "data/" << m_directoryStr << "/trajectoriesDir.dat";
        m_trajectoryDir.setPath(ss.str());
        if (m_append)
            m_trajectoryDir.openForAppend();
        else
            m_trajectoryDir.openForWrite();
        std::cout << "Writing trajectory orientations to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing trajectory orientations to file: " << ss.str() << std::endl;
    }
    
    if (constants::outputState)
    {
        ss.str("");
        ss << "data/" << m_directoryStr << "/recentState.dat";
        m_state.setPath(ss.str());
        std::cout << "Writing most recent state to file: " << ss.str() << std::endl;
        m_status.stream() << "Writing most recent state to file: " << ss.str() << std::endl;
    }
}

/* Reads simulation state from supplied path to a file
* - Assumes first line is simulation time
* - Assumes each subsequent line is particle data in the form of position, direction, next flip time
*/
bool io::readState(simulationState &a_state, const std::string a_path)
{
    file input;
    input.setPath(a_path);
    if(!input.openForRead())
        return false;
    std::string lineData;
    try
    {
        input.getline(lineData);
    }
    catch(const std::exception& e)
    {
        std::cerr << "Encountered exception: " << e.what() << std::endl;
        return false;
    }
    double clockTime;
    std::istringstream iss(lineData);
    if (!(iss >> clockTime))
    {
        std::cerr << "Error reading file " << a_path << std::endl;
        input.close();
        return false;
    }
    a_state.setSimulationClockTime(clockTime);
    double pos, timeGap;
    int32_t dir;
    uint32_t particlesRead = 0;
    std::vector<double> positions, flipTimes;
    std::vector<uint8_t> dirIdxs;
    for (uint32_t idx = 0; idx < constants::N; ++idx)
    {
        try
        {
            input.getline(lineData);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Encountered exception: " << e.what() << std::endl;
            break;
        }
        iss.clear();
        iss.str(lineData);
        if (!(iss >> pos >> dir >> timeGap))
        {
            break;
        }
        if ((pos < 0.0) || (pos > constants::L))
        {
            std::cerr << "Particle from file out of bounds." << std::endl;
            break;
        }
        positions.push_back(pos);
        uint16_t dirIdx = (dir > 0) ? static_cast<uint16_t>(1u) : static_cast<uint16_t>(0u);
        dirIdxs.push_back(static_cast<uint8_t>(dirIdx));
        flipTimes.push_back(timeGap);
        ++particlesRead;
    }
    input.close(); 
    if (particlesRead != constants::N)
    {   
        std::cerr << "Error reading file " << a_path << std::endl;
        return false;
    }
    a_state.setState(positions, dirIdxs, flipTimes);

    std::cout << "******************************************************" << std::endl;
    std::cout << "Read " << particlesRead << " particles from " << a_path << std::endl;
    std::cout << "Simulation restarting at clock time " << clockTime << std::endl;
    std::cout << "******************************************************" << std::endl << std::endl;

    m_status.stream() << "******************************************************" << std::endl;
    m_status.stream() << "Read " << particlesRead << " particles from " << a_path << std::endl;
    m_status.stream() << "Simulation restarting at clock time " << clockTime << std::endl;
    m_status.stream() << "******************************************************" << std::endl << std::endl;

    return true;    
}

/* 
* Returns true if this time step incurs a status output
*/
bool io::outputStatusCondition(const simulationState &a_state) const noexcept 
{
    return (!(a_state.simulationStep() % m_outputIntervalSteps));
}

/* 
* Returns true if this time step incurs a state output
*/
bool io::outputStateCondition(const simulationState &a_state) const noexcept 
{
    return (!(a_state.simulationStep() % m_outputStateIntervalSteps));
}

/* 
* Returns true if this time step incurs a data/field output
*/
bool io::outputDataCondition(const simulationState &a_state) const noexcept 
{
    return ((a_state.simulationTime() > constants::burn) && (!(a_state.simulationStep() % m_outputDataIntervalSteps)));
}

/* 
* Returns true if this time step incurs trajectory output
*/
bool io::outputTrajectoriesCondition(const simulationState &a_state) const noexcept 
{
    return ((a_state.simulationTime() > constants::burn) && (!(a_state.simulationStep() % m_outputTrajectoryIntervalSteps)));
}

/* 
* Output system status
* - checks if timestep is an output step
*/
void io::outputStatus(const simulationState &a_state) noexcept 
{
    if (outputStatusCondition(a_state))
    {
        auto tInterval = clock.interval();
        auto tDuration = clock.durationNow();
        m_status.stream() << "simulation run time: " << static_cast<double>(m_outputCount) * m_outputInterval
            << ", simulation clock time: " << a_state.simulationClockTime() 
            << ", duration since last status output (seconds): " << tInterval << ", total duration (seconds): " << tDuration << std::endl << std::flush;
        std::cout << "simulation time: " << static_cast<double>(m_outputCount) * m_outputInterval 
            << ", simulation clock time: " << a_state.simulationClockTime() 
            << ", duration since last status output (seconds): " << tInterval <<", total duration (seconds): "<< tDuration << std::endl << std::flush;
        
        /*
        const double averageSpeed  = a_state.averageSpeed();
        const double averageOrientation  = a_state.averageOrientation();

        m_status.stream() << "Average speed = " << averageSpeed << ", Average Orientation = " << averageOrientation << std::endl;
        std::cout << "Average speed = " << averageSpeed << ", Average Orientation = " << averageOrientation << std::endl;
        */
       
        ++m_outputCount;
    }
}

/* 
* Write simulation time to the supplied file
*/
void io::outputTime(const simulationState &a_state, file &a_file) noexcept 
{
    a_file.stream() << a_state.simulationClockTime() << std::endl;
}

/* 
* Write particle positions the supplied file
*/
void io::outputPositions(const simulationState &a_state, file &a_file) noexcept 
{
    for (std::size_t idx = 0u; idx + 1u < constants::N; ++idx)
        a_file.stream() << a_state.position(a_state.arrayPosOfParticle(idx)) << ",";
    if (constants::N > 0u) a_file.stream() << a_state.position(a_state.arrayPosOfParticle(constants::N - 1u)) << std::endl;
}

/* 
* Write particle directions the supplied file
*/
void io::outputDirections(const simulationState &a_state, file &a_file) noexcept 
{
    for (std::size_t idx = 0u; idx + 1u < constants::N; ++idx)
        a_file.stream() << a_state.velocity(a_state.arrayPosOfParticle(idx)) << ",";
    if(constants::N > 0u) a_file.stream() << a_state.velocity(a_state.arrayPosOfParticle(constants::N - 1u)) << std::endl;
}

/* 
* Write particle next tumble times the supplied file
*/
void io::outputFlipTimes(const simulationState &a_state, file &a_file) noexcept 
{
    for (std::size_t idx = 0u; idx + 1u < constants::N; ++idx)
        a_file.stream() << a_state.flipTime(a_state.arrayPosOfParticle(idx)) << ",";
    if(constants::N > 0u) a_file.stream() << a_state.flipTime(a_state.arrayPosOfParticle(constants::N - 1u)) << std::endl;
}

/* 
* Write particle state to the supplied file
*/
void io::outputParticleState(const simulationState &a_state, file &a_file, uint32_t a_idx) noexcept 
{
    a_file.stream() << std::setprecision(10u) << a_state.position(a_state.arrayPosOfParticle(a_idx)) << " "
        << a_state.velocity(a_state.arrayPosOfParticle(a_idx)) << " "
        << a_state.flipTime(a_state.arrayPosOfParticle(a_idx)) << std::endl;
}

/* 
* Output left and right moving particle density fields
*/
void io::outputSpecies(const std::vector<double> &a_rho, const std::vector<double> &a_psi) noexcept 
{
    const std::size_t rhoLen = a_rho.size();
    for (std::size_t idx = 0u; idx + 1u < rhoLen; ++idx)
        m_leftField.stream() << std::setprecision(10u) << 0.5 * (a_rho[idx] - a_psi[idx]) << ",\t";
    if (rhoLen) m_leftField.stream() << 0.5 * (a_rho[rhoLen - 1u] - a_psi[rhoLen - 1u]) << std::endl;              
    for (std::size_t idx = 0u; idx + 1u < rhoLen; ++idx)
        m_rightField.stream() << std::setprecision(10u) << 0.5 * (a_rho[idx] + a_psi[idx]) << ",\t";
    if (rhoLen) m_rightField.stream() << 0.5 * (a_rho[rhoLen - 1u] + a_psi[rhoLen - 1u]) << std::endl;              
}

/* 
* Output particle density field
*/
void io::outputField(const std::vector<double> &a_rho) noexcept 
{
    const std::size_t rhoLen = a_rho.size();
    for (std::size_t idx = 0u; idx + 1u < rhoLen; ++idx)
        m_densityField.stream() << std::setprecision(10) << a_rho[idx] << ",\t";
    if(rhoLen) m_densityField.stream() << a_rho[rhoLen-1] << std::endl;
}

/* 
* Output particle polarisation density field
*/
void io::outputPolarisationField(const std::vector<double> &a_psi) noexcept 
{
    const std::size_t psiLen = a_psi.size();
    for (std::size_t idx = 0u; idx + 1u < psiLen; ++idx)
        m_dirField.stream() << std::setprecision(10u) << a_psi[idx] << ",\t";
    if (psiLen) m_dirField.stream() << a_psi[psiLen - 1u] << std::endl;              
}

/* 
* Output particle density fields
* - includes logic to output different fields dependent on parameters
*/
void io::outputFields(const simulationState &a_state) noexcept 
{
    const double &fieldScale = m_fieldScale; 

    if (constants::outputSpeciesField)
    {
        const std::vector<double> &rho = a_state.densityField(fieldScale);
        const std::vector<double> &psi = a_state.polarisationDensityField(fieldScale);
        outputSpecies(rho,psi);
        if (constants::outputField)        
            outputField(rho);
        if (constants::outputDirField)
            outputPolarisationField(psi);
    }
    else
    {
        if (constants::outputField)        
        {
            const std::vector<double> &rho = a_state.densityField(fieldScale);
            outputField(rho);
        }
        if (constants::outputDirField)
        {
            const std::vector<double> &psi = a_state.polarisationDensityField(fieldScale);
            outputPolarisationField(psi);
        }
    }
}

/* 
* Output particle trajectories
*/
void io::outputTrajectories(const simulationState &a_state) noexcept 
{
    outputTime(a_state, m_trajectoryTime);
    outputPositions(a_state, m_trajectoryPos);
    outputDirections(a_state, m_trajectoryDir);
}

/* 
* Output simulation state
*/
void io::outputState(const simulationState &a_state) noexcept 
{
    m_state.openForWrite();
    outputTime(a_state, m_state);
    for (uint32_t idx = 0; idx < constants::N; ++idx)
        outputParticleState(a_state, m_state, idx);
    m_state.close();
}

/* 
* Output data
*/
void io::outputData(const simulationState &a_state) noexcept 
{
    if (constants::outputField || constants::outputDirField || constants::outputSpeciesField)
    {
        if (outputDataCondition(a_state))
            outputFields(a_state);        
    }

    if (constants::outputState)
    {
        if (outputStateCondition(a_state))
            outputState(a_state);
    }

    if (constants::outputTrajectories)
    {
        if (outputTrajectoriesCondition(a_state))
            outputTrajectories(a_state);
    }
}

/* 
* Write parameter values to file
*/
void io::outputParameters() noexcept 
{
    std::cout << "writing parameters to file: " << m_parameter.path() << std::endl << std::endl;
    m_status.stream() << "writing parameters to file: " << m_parameter.path() << std::endl << std::endl;

    m_parameter.openForWrite();
    m_parameter.stream() << "Simulation started at: " << m_timestampStr << std::endl;
    m_parameter.stream() << "Data written to directory: " << m_fullDirectoryStr << std::endl << std::endl;
    
    m_parameter.stream() << "SIMULATION PARAMETERS:" << std::endl << std::endl;
    
    m_parameter.stream() << "N [constants::N] = " << constants::N << std::endl;
    m_parameter.stream() << "L [constants::L] = " << constants::L << std::endl;
    m_parameter.stream() << "dt (seconds) [constants::dt] = " << constants::dt << ", (const sqrt(dt) = " << constants::sqrtDt << ")" << std::endl;
    m_parameter.stream() << "interaction distance [constants::interactDist] = " << constants::interactDist << std::endl;
    m_parameter.stream() << "interval between data sorts (time steps) [constants::sortInterval] = " << constants::sortInterval << " (every" << constants::sortIntervalSteps << " integration steps)" <<std::endl;
    m_parameter.stream() << "number of threads [constants::numThreads] = " << constants::numThreads << std::endl << std::endl;

    m_parameter.stream() << "OUTPUT OPTIONS/PARAMETERS:" << std::endl << std::endl;

    m_parameter.stream() << "Output density field [constants::outputField]: " << ((constants::outputField) ? "yes/true" : "no/false") <<std::endl;
    m_parameter.stream() << "Output polarisation density field [constants::outputDirField]: " << ((constants::outputDirField) ? "yes/true" : "no/false") <<std::endl;
    m_parameter.stream() << "Output species (left/right) density fields [constants::outputSpeciesField]: " << ((constants::outputSpeciesField) ? "yes/true" : "no/false") <<std::endl;
    m_parameter.stream() << "Output particle trajectories [constants::outputTrajectories]: " << ((constants::outputTrajectories) ? "yes/true" : "no/false") <<std::endl;
    m_parameter.stream() << "field output resolution [constants::fieldResolution] = " << constants::fieldResolution << std::endl;
    m_parameter.stream() << "field basis-function scale/width [constants::fieldScale] = " << m_fieldScale << std::endl << std::endl;
    
    m_parameter.stream() << "DURATION/OUTPUT INTERVAL PARAMETERS" << std::endl << std::endl;

    m_parameter.stream() << "burn time (seconds) [constants::burn] = " << constants::burn << std::endl;
    m_parameter.stream() << "total simulation time (seconds) [constants::totalTime] = " << constants::totalTime << std::endl;
    m_parameter.stream() << "progress output interval (seconds) [constants::outputInterval] = " << m_outputInterval << std::endl;
    m_parameter.stream() << "data output interval (seconds) [constants::outputDataInterval] = " << m_outputDataInterval << std::endl;
    m_parameter.stream() << "state output interval (seconds) [constants::outputStateInterval] = " << m_outputStateInterval << std::endl;
    m_parameter.stream() << "trajectory output interval (seconds) [constants::outputTrajInterval] = " << m_outputTrajectoryInterval << std::endl << std::endl;

    m_parameter.stream() << "SDE/DENSITY/QUORUM DEFINITIONS/CONVENTIONS:" << std::endl << std::endl;
    m_parameter.stream() << "model density 'constants::rho' = one particle probability density / constants::rhoC" << std::endl;
    m_parameter.stream() << "constants::rhoC = prob. density that arrests particle propulsion" << std::endl;
    m_parameter.stream() << "raw one particle probability density = 1 / L = 'RHO' * 'RHO_C'" << std::endl;
    m_parameter.stream() << "physical density (N / L) = N * one particle probability density" << std::endl;
    m_parameter.stream() << "total quorum sensing effect = quorum strengh * local physical density = local physical density / (rhoC * N) = local prob density / rhoC" << std::endl;
    m_parameter.stream() << "SDE: dx = propulsion_strength * (1 - quorum_strength * local_density) * dt + noise_strength * dW" << std::endl;
    m_parameter.stream() << "local_density = local physical density : N_loc / V_loc; V_loc = 2.0 * constants::interactDist" << std::endl << std::endl;

    const double diffusionStrength = 0.5 * (constants::noiseStrength) * (constants::noiseStrength);
    const double gammaBar = 2.0 * (constants::gammaMinus) * (constants::gammaPlus) / (constants::gammaMinus + constants::gammaPlus);

    m_parameter.stream() << "MODEL PARAMETERS - SPECIFIED:" << std::endl << std::endl;
    
    m_parameter.stream() << "uniform state model density [constants::rho] = " << constants::rho << std::endl;
    m_parameter.stream() << "propulsion strength [constants::propulsionStrength] = " << constants::propulsionStrength << std::endl;
    m_parameter.stream() << "noise strength (sqrt(2D)) [constants::noiseStrength] = " << constants::noiseStrength << std::endl;
    m_parameter.stream() << "beta (symmetry parameter) [constants::beta] = " << constants::beta << std::endl;
    m_parameter.stream() << "Peclet number (propulsion_strength / sqrt(mean_tumbling_rate * D)) [constants::Pe] = " << constants::Pe << std::endl << std::endl;

    m_parameter.stream() << "MODEL PARAMETERS - CALCULATED/AGGREGATE:" << std::endl << std::endl;

    m_parameter.stream() << "quorum sensing density constant constants::rhoC = (1 / (constants::rho * constants::L)) [constants::rhoC] = " << constants::rhoC << std::endl;
    m_parameter.stream() << "quorum sensing strength (1/(N * RHO_C)) [constants::quorumStrength] = " << constants::quorumStrength << std::endl;
    m_parameter.stream() << "diffusion strength, D = " << diffusionStrength << std::endl;
    m_parameter.stream() << "bulk tumbling rate, gamma_bar = (2 * gamma_plus * gamma_minus / (gamma_plus + gamma_minus)) = " << gammaBar << std::endl;
    m_parameter.stream() << "gamma_plus (rate of -ve -> +ve tumbles) [constants::gammaPlus] = " << constants::gammaPlus << std::endl;
    m_parameter.stream() << "gamma_minus (rate of +ve -> -ve tumbles) [constants::gammaMinus] = " << constants::gammaMinus << std::endl << std::endl;
    m_parameter.stream() << "Recalculated Pe number (propulsion_strength / sqrt(gamma_bar * D)) = " << constants::propulsionStrength / sqrt( diffusionStrength * gammaBar ) << std::endl;
    m_parameter.stream() << "Recalculated symmetry, beta = gamma_minus / gamma_plus = " << (constants::gammaMinus) / (constants::gammaPlus) << std::endl << std::endl;    

    m_parameter.stream() << "INITIAL CONDITIONS:" << std::endl << std::endl;

    if (constants::nucleate)
    {
        m_parameter.stream() << "Intial condition: top hat density function." << std::endl;
        m_parameter.stream() << "average model density [constants::rhoMips (= constants::rho)] = " << constants::rho << std::endl;
        m_parameter.stream() << "average equivalent probability density [constants::rhoMips * constants::rhoC] = " << constants::rho * constants::rhoC << std::endl;
        m_parameter.stream() << "initialised with model densitities: rho_H [constants::rhoMipsHigh] = " << constants::rhoMipsHigh << ", rho_L [constants::rhoMipsLow] = " << constants::rhoMipsLow << std::endl;
        m_parameter.stream() << "equivalent one particle probability densities (* constants::rhoC): " << (constants::rhoMipsHigh) * (constants::rhoC) << ", " << (constants::rhoMipsLow) * (constants::rhoC) << std::endl;
    }
    else{
        m_parameter.stream() << "Intial condition: constant density profile." << std::endl;
        m_parameter.stream() << "constant model density [constants::rho] = " << constants::rho << std::endl;
        m_parameter.stream() << "equivalent one particle probability density [constants::rho * constants::rhoC]: " << constants::rho * constants::rhoC << std::endl;
    }

    m_parameter.close();   
}
