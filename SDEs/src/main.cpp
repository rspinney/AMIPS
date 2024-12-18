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

#include <string>
#include "model.hpp" 

/*
    Up to 2 arguments:
    0 Arguments: will create a directory ./data/default/ and output files to it
    1 Argument:  will create a directory ./data/$1/ and output files to it
    2 Arguments: will create a directory ./data/$1/ and output files to it
                 will read in initial particle state from file at $2

*/

int main(int argc, char *argv[]) 
{
    const std::string folderStr = (argc > 1) ? std::string(argv[1u]) : "default";
    const std::string storedStateStr = (argc > 2) ? std::string(argv[2u]) : "";
    model simulation(folderStr, storedStateStr);
    simulation.run();
    return 0;
}