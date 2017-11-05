set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5.1)
# set( CMAKE_VERBOSE_MAKEFILE on )
set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_CXX_COMPILER clang )
set(THREADS_PREFER_PTHREAD_FLAG ON)

project(ImmunoDBSCAN)
# The version number.
set (ImmunoDBSCAN_VERSION_MAJOR 1)
set (ImmunoDBSCAN_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/ImmunoDBSCAN.h.in"
  "${PROJECT_BINARY_DIR}/cpp/ImmunoDBSCAN.h"
  )
  
#************************************************************  
# This variable has to be set manually!!!!  
set (CUDA_TOOLKIT_ROOT_DIR  /usr/local/cuda-8.0)

#************************************************************

find_package(Boost REQUIRED system program_options)
find_package(CUDA REQUIRED)
find_package(Threads REQUIRED)

include_directories(${Boost_INCLUDE_DIRS}, ${CUDA_INCLUDE_DIRS})
CUDA_INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}, ${CUDA_INCLUDE_DIRS})
#SET(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-Wall -Werror -pthread -std=c++11 -Wl,--no-as-needed" -D_MWAITXINTRIN_H_INCLUDED)
#SET(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS} "-Wall -Werror -pthread -std=c++11 -Wl,--no-as-needed")

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -Wextra -g -O2")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -O2")

SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_50,code=sm_52 --std=c++11)

#cuda_compile(BITSETVECTOR_O cuda/BitsetVector.cu)
cuda_compile(INTERLINK_O cuda/Interlink.cu)
cuda_add_executable(psdbscan cpp/psdbscan.cpp ${INTERLINK_O})
#add_executable(main cpp/psdbscan.cpp ${InterlinkClass_O})
#target_link_libraries(main ${Boost_LIBRARIES} ${CUDA_LIBRARIES} )

target_link_libraries(psdbscan Threads::Threads ${Boost_LIBRARIES} ${CUDA_LIBRARIES})
#target_link_libraries(main Threads::Threads ${Boost_LIBRARIES})

# cuda_compile(CORETEST_O cuda/IsCorePoint.cu)
# cuda_add_executable(coretest ${CORETEST_O} )
# target_link_libraries(coretest ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY})
