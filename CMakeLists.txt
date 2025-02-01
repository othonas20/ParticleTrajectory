cmake_minimum_required(VERSION 3.23.2)
project(ParticleTrajectory)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# include Eigen
set(Eigen_dir
        /home/vs/eigen-3.4.0
#         C:/Users/VS0121/eigen-3.4.0
        )
include_directories(${Eigen_dir})

# Include subdirectories for libraries
include(${CMAKE_CURRENT_LIST_DIR}/lib/CMakeLists.txt)

# Include header files
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)

# include source files
include(src/CMakeLists.txt)

# Add the executable for the main application
add_executable(
        ParticleTrajectory ${LIB_SRC} ${SRC}
        src/exerciseA.cpp
        )
