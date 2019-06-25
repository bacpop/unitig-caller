cmake_minimum_required(VERSION 3.4)
project(map_strings)
set(CMAKE_CXX_STANDARD 17)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/seqan3/build_system")
find_package (SeqAn3 REQUIRED)

find_package(pybind11 REQUIRED)
pybind11_add_module(map_strings src/map_strings.cpp src/map_bindings.cpp)

target_link_libraries(map_strings PRIVATE seqan3::seqan3 pthread)
