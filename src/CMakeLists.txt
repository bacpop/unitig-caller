cmake_minimum_required(VERSION 3.0)
project(map_strings)

find_package(pybind11 REQUIRED)
pybind11_add_module(map_strings fasta.cpp map_strings.cpp map_bindings.cpp)
