set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5.1)

project(ImmunoDBSCAN)
# The version number.
set (ImmunoDBSCAN_VERSION_MAJOR 1)
set (ImmunoDBSCAN_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/ImmunoDBSCAN.h.in"
  "${PROJECT_BINARY_DIR}/cuda/ImmunoDBSCAN.h"
  )

find_package(Boost REQUIRED system program_options)
find_package(CUDA REQUIRED)

include_directories(${Boost_INCLUDE_DIRS}, ${CUDA_INCLUDE_DIRS})

add_executable(main cuda/HelloWorld.cpp)
target_link_libraries(main ${Boost_LIBRARIES} ${CUDA_LIBRARIES})

SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_50,code=sm_52 --std=c++11)
CUDA_INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}, ${CUDA_INCLUDE_DIRS})
cuda_compile(CORETEST_O cuda/IsCorePoint.cu)
cuda_add_executable(coretest ${CORETEST_O} )
target_link_libraries(coretest ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY})
