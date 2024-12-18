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

#ifndef CONSTFUNCS_H
#define CONSTFUNCS_H

#include <limits>   

namespace Implementation
{
double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
{
    return curr == prev
        ? curr
        : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}
}

/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
*/
double constexpr sqrtConst(double x)
{
    return (x >= 0 && x < std::numeric_limits<double>::infinity())
        ? Implementation::sqrtNewtonRaphson(x, x, 0)
        : std::numeric_limits<double>::quiet_NaN();
}

#endif /*CONSTFUNCS*/