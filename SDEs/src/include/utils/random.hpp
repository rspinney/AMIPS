/*****************************************************************************/
/***************** Copyright (C) 2020-2021, Richard Spinney. *****************/
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

#ifndef RANDOM_H
#define RANDOM_H

/*
    class wrapper for c++11 uniform distribution of psuedo-random number
	using templates for different generators, e.g. Mersenne Twister algo.

	possible types for template argument:
	std::minstd_rand
	std::mt19937
	std::mt19937_64

	default template argument requires c++17
*/

#include <random>     // std::random_device, std::normal_distribution etc
#include <array>      // std::array
#include <cstdint>    // uint32_t etc
#include <algorithm>  // std::generate_n

template<typename T = std::mt19937> //default type - >= c++17 only
class randNumGen
{
	T m_eng;
	std::uniform_real_distribution<double> m_uniform;
	std::normal_distribution<double> m_gaussian;
	void init_seed()
	{
		std::array<int32_t, 624u> seedData;//fill the MT registers with random numbers
		std::random_device rndDev;
		std::generate_n(seedData.data(), seedData.size(), std::ref(rndDev));
		std::seed_seq seq(std::begin(seedData), std::end(seedData));
		m_eng = T(seq);
	}
public:
	randNumGen() : m_uniform(0.0L,1.0L), m_gaussian(0.0L,1.0L){init_seed();}
	randNumGen(const uint32_t a_seed): m_eng(a_seed), m_uniform(0.0L, 1.0L), m_gaussian(0.0L, 1.0L){};
	double uniform(){return m_uniform(m_eng);}
	double gaussian(){return m_gaussian(m_eng);}
};

#endif /*RANDOM_H*/
