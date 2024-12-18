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

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <array>    // std::array
#include <vector>   // std::vector

#include "random.hpp"

// #ifdef __cpp_lib_hardware_interference_size
//     using std::hardware_constructive_interference_size;
//     using std::hardware_destructive_interference_size;
// #else
    //64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
    constexpr std::size_t hardware_constructive_interference_size = 256;
    constexpr std::size_t hardware_destructive_interference_size = 256;
//#endif

struct integrationInfo
{
    uint32_t particle = 0u, threadId = 0u, bin = 0u;
};

struct particleAndBinPair
{
    uint32_t particle, binOld, binNew;
    particleAndBinPair() : particle(0u), binOld(0u), binNew(0u){}
    particleAndBinPair(const uint32_t a_particle, const uint32_t a_binOld, const uint32_t a_binNew) :
        particle(a_particle), binOld(a_binOld), binNew(a_binNew){}
};

class alignas(hardware_destructive_interference_size) 
threadObjects
{
    using rngGen = std::mt19937_64; 
    randNumGen<rngGen> m_rng;
    std::array<uint32_t,2u> m_setBounds;
    std::vector<particleAndBinPair> m_updateList;
    public:
    threadObjects(const uint32_t a_threadIdx) : m_rng(a_threadIdx + 1234u), m_setBounds({0,0}){}
    threadObjects() : m_rng(){}
    randNumGen<rngGen>& rng(){ return m_rng;}
    std::array<uint32_t,2u>& setBounds(){ return m_setBounds;}
    std::vector<particleAndBinPair>& updateList(){ return m_updateList;}
};

#endif /*STRUCTURES_H*/
