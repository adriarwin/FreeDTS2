cmake_minimum_required(VERSION 3.10)

# Project name and version
project(FreeDTS2 VERSION 1.0)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define source directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/dts_src)
set(CONVERT_DIR ${CMAKE_SOURCE_DIR}/dts_convert)
set(GENERATE_DIR ${CMAKE_SOURCE_DIR}/dts_generate)

# Gather all source files
file(GLOB SRC_FILES ${SRC_DIR}/*.cpp)
file(GLOB CONVERT_FILES ${CONVERT_DIR}/*.cpp)
file(GLOB GENERATE_FILES ${GENERATE_DIR}/*.cpp)

# Include directories (add if needed)
include_directories(${SRC_DIR} ${CONVERT_DIR} ${GENERATE_DIR})

# Define options for building with OpenMP or MPI
option(USE_OPENMP "Use OpenMP" OFF)
option(USE_MPI "Use MPI" OFF)

# Create targets for different build types
if(USE_MPI)
    find_package(MPI REQUIRED)
    add_executable(DTS ${SRC_FILES})
    target_link_libraries(DTS MPI::MPI_CXX)
    target_compile_definitions(DTS PRIVATE MPI_DETECTED)
else()
    add_executable(DTS ${SRC_FILES})
endif()

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(DTS OpenMP::OpenMP_CXX)
endif()

add_executable(CNV ${CONVERT_FILES})
add_executable(GEN ${GENERATE_FILES})

# Install targets to the bin directory
install(TARGETS DTS CNV GEN DESTINATION bin)
