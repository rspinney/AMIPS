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
    Functions to get permutation of indices corresponding to sorted vector,
    And to apply that permuation in place to given vectors of the same length
*/

#ifndef SORTPERM_H
#define SORTPERM_H

#include <vector>    // std::vector
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::swap
//#include <execution> // std::execution::par_unseq

template <typename T, typename Compare>
inline std::vector<std::size_t> sortedPermutation(const std::vector<T>& a_vec, const Compare& a_compare = std::less<T>())
{
    std::vector<std::size_t> permutation(a_vec.size());
    std::iota(permutation.begin(), permutation.end(), 0u);
    // std::sort(std::execution::par_unseq, permutation.begin(), permutation.end(), // valgrind says this leaks...
    //     [&](std::size_t i, std::size_t j){ return a_compare(a_vec[i], a_vec[j]); });
    std::sort(permutation.begin(), permutation.end(),
        [&](std::size_t i, std::size_t j){ return a_compare(a_vec[i], a_vec[j]); });
    return permutation;
}

template <typename T>
inline void applyPermutation(std::vector<T>& a_vec, const std::vector<std::size_t>& a_permutation)
{
    std::vector<bool> done(a_vec.size());
    for (std::size_t i = 0u; i < a_vec.size(); ++i)
    {
        if (done[i])
            continue;
        done[i] = true;
        std::size_t jPrev = i;
        std::size_t j = a_permutation[i];
        while (i != j)
        {
            std::swap(a_vec[jPrev], a_vec[j]);
            done[j] = true;
            jPrev = j;
            j = a_permutation[j];
        }
    }
}

#endif /* SORTPERM_H */
