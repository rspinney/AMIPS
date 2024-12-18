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
    class to encapsulate conversion of individual particle
    positions into a continuous field.

    uses gaussian basis functions at given parametrised scale

    basis functions are (must be) pre-computed using gaussianBasis::set(scale)

    parallelism achieved by the interface provided by "parallel.hpp"
*/

#ifndef FIELD_H
#define FIELD_H

#include <vector>       // std::vector
#include <algorithm>    // std::fill
#include <cmath>        // std::exp, std::fmod
#include <cstdint>      // uint32_t etc.

#include "parameters.hpp"
#include "staticasserts.hpp"
#include "parallel.hpp"

class field
{
    const double m_resolution;         // size of discretised bin
    const double m_length;             // size of field in units L
    const std::size_t m_numBins;       // nummber of discretised elements of field
    std::vector<double> m_field;       // field data buffer
    double m_scaling;                  // multiplier for output to adjust from the native probability density
    
    // threading objects
    const uint32_t m_numThreads;
    parallel m_parallel; // parallelise calc
    std::vector<std::vector<double> > m_fieldThread; //independent buffers for work threads, false sharing unlikely
    
    class gaussianBasis //save lengthscale state for Gaussian basis functions
    {
        std::vector<double> m_cache; // precompute basis function at range of distances
        double m_prefact, m_expConst; // gaussian parameters
        int64_t m_offset; // index of most extremal position considered from particle
        const double m_maxDistMult, m_resolution;
        
        double value(const double &a_dx) const noexcept 
        {
            return m_prefact * std::exp(- m_expConst * a_dx * a_dx);
        }

        void cache(double a_scale) noexcept 
        {
            m_cache.clear();
            m_offset = static_cast<int64_t>(m_maxDistMult * a_scale / m_resolution);
            m_cache.reserve(2u * m_offset);
            for (int64_t i = -m_offset; i < m_offset; ++i)
            {
                const double dx = static_cast<double>(i) * m_resolution;
                m_cache.emplace_back(value(dx));                
            }
        }
        
        std::size_t PBCwrap(const int64_t a_idx, const int64_t a_maxIdx) const noexcept 
        {
            return ((a_idx % a_maxIdx) + a_maxIdx) % a_maxIdx;
        }

        public:
        
        gaussianBasis(const double a_maxDistMult, const double a_resolution) :
            m_maxDistMult(a_maxDistMult),
            m_resolution(a_resolution)
        {
            set(1.0L);
        }
        
        void set(double a_sigma) noexcept 
        {
            m_prefact = 1.0 / std::sqrt(2.0 * a_sigma * a_sigma * constants::pi);
            m_expConst = 0.5 / (a_sigma * a_sigma);
            cache(a_sigma);
        }

        void addBasis(std::vector<double> & a_buffer, double a_x, double a_v = 1.0L) const noexcept 
        {
            const int64_t xIndex = static_cast<int64_t>(a_x / m_resolution);
            for (int64_t i = 0L; i < 2L * m_offset; ++i)
            {
                const std::size_t fieldIndex = PBCwrap(i + xIndex - m_offset, static_cast<int64_t>(a_buffer.size()));
                a_buffer[fieldIndex] +=  a_v * m_cache[i];
            }
        }
    };

    gaussianBasis m_gaussian;

    public:

    field(const uint32_t a_numThreads, double a_resolution, double a_length, double a_scaling) :
        m_resolution(a_resolution),
        m_length(a_length),
        m_numBins(static_cast<size_t>(m_length / m_resolution)),
        m_field(m_numBins, 0.0L),
        m_scaling(a_scaling),
        m_numThreads(a_numThreads),
        m_parallel(a_numThreads),
        m_fieldThread(m_numThreads, std::vector<double>(m_numBins,0.0L)),
        m_gaussian(5.0L, m_resolution) // include gaussian up to first argument (e.g. 5.0L) times sigma from the mean
    {}
    
    // field calculations -
    // outputs a field formed of Gaussian basis functions
    // for simple density field \sum_i \rho_i * dx = num_particles
    // i.e. outputs density field, not probability or particle fraction field (or whatever it might be called)

    // shared function to gather results from all threads

    void sumFieldThread(const uint32_t a_threadIdx, const std::size_t a_N) noexcept 
    {
        const uint32_t idxPerThread = (static_cast<uint32_t>(m_numBins) / m_numThreads);
        const uint32_t start = a_threadIdx * idxPerThread;
        const uint32_t end = (a_threadIdx + 1u < m_numThreads) * ((a_threadIdx + 1u) * idxPerThread) + (a_threadIdx + 1u >= m_numThreads) * (static_cast<uint32_t>(m_numBins));
        const double invN = 1.0 / static_cast<double>(a_N);
        for (uint32_t idx = start; idx < end; ++idx)
        {
            double sum = 0.0L;
            for (uint16_t threadIdx = 0u; threadIdx < m_numThreads; ++threadIdx)
                sum += m_fieldThread[threadIdx][idx];
            m_field[idx] = m_scaling * invN * sum;
        }
    }

    // density field

    void calcFieldThreaded(const uint32_t a_threadIdx, const std::vector<double> &a_data) noexcept 
    {
        std::vector<double> &fieldThread = m_fieldThread[a_threadIdx];
        std::fill(fieldThread.begin(), fieldThread.end(), 0.0L);
        const uint32_t dataSize = static_cast<uint32_t>(a_data.size());
        const uint32_t idxPerThread = (dataSize / m_numThreads);
        const uint32_t start = a_threadIdx * idxPerThread;
        const uint32_t end = (a_threadIdx + 1u < m_numThreads) * ((a_threadIdx + 1u) * idxPerThread) + (a_threadIdx + 1u >= m_numThreads) * (dataSize);
        for (uint32_t idx = start; idx < end; ++idx)
            m_gaussian.addBasis(fieldThread, a_data[idx]);
        m_parallel.synchronise();
        sumFieldThread(a_threadIdx, a_data.size());
    }

    const std::vector<double>& calcField(const std::vector<double> &a_data, const double &a_scale) noexcept 
    {
        m_gaussian.set(a_scale);
        m_parallel.parallelIndexedTask([&](uint32_t a_threadIdx){calcFieldThreaded(a_threadIdx,std::cref(a_data));});
        return m_field;
    }

    //polarisation density field

    void calcVFieldThreaded(uint32_t a_threadIdx, const std::vector<double> &a_data, const std::vector<uint8_t> &a_vIndex, const std::array<double,2> &a_vVals) noexcept 
    {
        std::vector<double> &fieldThread = m_fieldThread[a_threadIdx];
        std::fill(fieldThread.begin(), fieldThread.end(), 0.0L);
        const uint32_t dataSize = static_cast<uint32_t>(a_data.size());
        const uint32_t idxPerThread = (dataSize / m_numThreads);
        const uint32_t start = a_threadIdx * idxPerThread;
        const uint32_t end = (a_threadIdx + 1u < m_numThreads) * ((a_threadIdx + 1u) * idxPerThread) + (a_threadIdx + 1u >= m_numThreads) * (dataSize);
        for (uint32_t idx = start; idx < end; ++idx)
            m_gaussian.addBasis(fieldThread, a_data[idx], a_vVals[a_vIndex[idx]]);
        m_parallel.synchronise();
        sumFieldThread(a_threadIdx, a_data.size());
    }

    const std::vector<double>& calcField(const std::vector<double> &a_data, const std::vector<uint8_t> &a_vIndex, const std::array<double,2> &a_vVals, const double &a_scale) noexcept 
    {
        m_gaussian.set(a_scale);
        m_parallel.parallelIndexedTask([&](uint32_t a_threadIdx){calcVFieldThreaded(a_threadIdx,std::cref(a_data),std::cref(a_vIndex),std::cref(a_vVals));});
        return m_field;
    }
};

#endif /* FIELD_H */
