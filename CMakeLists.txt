cmake_minimum_required(VERSION 3.4)
set(TARGET_NAME unitig_query)
project("${TARGET_NAME}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

IF(CMAKE_COMPILER_IS_GNUCC OR "$ENV{SKETCHLIB_INSTALL}" STREQUAL "conda")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF()

# Vendored libraries
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vendor)

#Find include dirs for z and bifrost
FIND_PATH(z_INCLUDE_DIR NAMES zlib.h HINTS ${CONDA_DIR}/include)
FIND_PATH(bifrost_INCLUDE_DIR bifrost/CompactedDBG.hpp HINTS ${CONDA_DIR}/include)

# Find libraries for z and bifrost
FIND_LIBRARY(z_LIBRARY NAMES z HINTS ${CONDA_DIR}/lib REQUIRED)
FIND_LIBRARY(bifrost_LIBRARY NAMES bifrost HINTS ${CONDA_DIR}/lib REQUIRED)

# Set package found as true
FIND_PACKAGE_HANDLE_STANDARD_ARGS(z DEFAULT_MSG z_INCLUDE_DIR z_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(bifrost DEFAULT_MSG bifrost_INCLUDE_DIR bifrost_LIBRARY)

# Set dirs if multiple include/libraries found
SET( z_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( z_LIBRARIES ${z_LIBRARY} )
SET( bifrost_INCLUDE_DIRS ${z_INCLUDE_DIR} )
SET( bifrost_LIBRARIES ${bifrost_LIBRARY} )

execute_process(COMMAND pybind11-config --cmakedir OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_DIR)
find_package(pybind11 2.6 CONFIG REQUIRED)
target_sources("${TARGET_NAME}" src/map_strings.cpp src/map_bindings.cpp src/bifrost.cpp)

if(UNIX AND (NOT APPLE OR NOT CMAKE_COMPILER_IS_GNUCC))
    set_target_properties("${TARGET_NAME}" PROPERTIES
                          INTERPROCEDURAL_OPTIMIZATION ON)
endif()

target_link_libraries("${TARGET_NAME}" PRIVATE pybind11::module ${bifrost_LIBRARIES} ${z_LIBRARIES})
if(DEFINED ENV{CONDA_PREFIX} AND (NOT APPLE OR CMAKE_COMPILER_IS_GNUCC OR ENV{UC_INSTALL} EQUAL "conda"))
    target_link_libraries("${TARGET_NAME}" PRIVATE gomp)
else()
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries("${TARGET_NAME}" PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()
