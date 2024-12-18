cmake_minimum_required(VERSION 3.10)

# get rid of invasive conda path nonsense

if(DEFINED ENV{CONDA_PREFIX})
	message(VARR="${CONDA_PREFIX}")
  set(h5_ignore_path
    $ENV{CONDA_PREFIX}/bin $ENV{CONDA_PREFIX}/lib $ENV{CONDA_PREFIX}/include
    $ENV{CONDA_PREFIX}/Library/bin $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/Library/include
  )
  list(APPEND CMAKE_IGNORE_PATH ${h5_ignore_path})
endif()

# set the project name
set(EXEC_NAME AMIPS_SDEs)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(${EXEC_NAME}  VERSION 0.1 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#compiler flags - tested on OSX and Linux
if (MSVC)
	SET(CMAKE_CXX_FLAGS_RELEASE "/O2")
	SET(CMAKE_C_FLAGS_RELEASE "/O2")
elseif(APPLE) #assuming clang
	set(CMAKE_CXX_FLAGS_DEBUG "-g3 -flto -ggdb -Wfatal-errors")
	set(CMAKE_CXX_FLAGS_RELEASE "-flto=auto -Ofast -Wall -Wextra -Wpedantic -ftree-vectorize -march=native")
	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
else() #assuing gcc
	set(CMAKE_CXX_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address -g3 -ggdb -O0 -Wfatal-errors") #-Weffc++ -Werror
	set(CMAKE_CXX_FLAGS_RELEASE "-flto=auto -Ofast -fuse-linker-plugin -Werror -Wfatal-errors  -Wall -Wextra -Wpedantic -Wcast-align -Wunused -Wold-style-cast -Wpointer-arith -Wcast-qual -Wno-missing-braces -Wconversion -Wshadow -ftree-vectorize -finline-functions -funroll-loops -msse2 -mfpmath=sse -march=native")
	SET(CMAKE_EXE_LINKER_FLAGS "-flto=auto -Ofast -fuse-linker-plugin")
	SET(CMAKE_SHARED_LINKER_FLAGS "-flto=auto -Ofast -fuse-linker-plugin")
	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-flto=auto -Ofast -fuse-linker-plugin")	
endif()
    
set(MAIN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# sources
set(SRCS ./src/main.cpp 
		./src/model.cpp 
		./src/integrator.cpp 
		./src/io.cpp 
		./src/simulationstate.cpp
		./src/binstructure.cpp)

# add the executable
add_executable(${EXEC_NAME} ${SRCS})

#headers
target_include_directories(${EXEC_NAME}  PRIVATE 
		./src/include/AMIPS 
		./src/include/utils
		)


find_package(Threads REQUIRED)
target_link_libraries(${EXEC_NAME} PUBLIC Threads::Threads)

