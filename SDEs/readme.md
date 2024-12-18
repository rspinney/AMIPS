# Many body simulation of SDEs representing biased QSAPs

## Installation

You can generate the build files and build the executable using the supplied CMake script using

`mkdir build && cd build && cmake .. && make`

This will build an executable file `AMIPS_SDEs` in the `build` directory.

## Parameters 

To change parameters and initial conditions edit the values in `./src/include/AMIPS/parameters.hpp` and recompile with `make`.

## Output

Depending on parameter values, various coarse grained empirical fields are output every `outputDataInterval` simulation seconds to files in the specified data path. If `outputState == true` then the total system state will be periodically written out to `recentState.dat` in the data path every `outputStateInterval`.

## Running

Running the executable with one argument e.g.

`./AMIPS_SDEs directory`

will start a simulation with output directed to files in

`./data/directory/`

Running the executable with two arguments

`./AMIPS_SDEs directory data/directory/recentState.dat`

will start a simulation with output directed to files in

`./data/directory/`

and initialise the particles according to data found in

`./data/directory/recentState.dat`

appending further output to any data files already in the data directory.