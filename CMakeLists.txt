cmake_minimum_required(VERSION 3.0)
project(map_strings)

find_package(pybind11 REQUIRED)
pybind11_add_module(map_strings src/fasta.cpp src/map_strings.cpp src/map_bindings.cpp)
target_link_libraries(map_strings PRIVATE pthread)
