cmake_minimum_required(VERSION 3.0)
project(map_strings)

find_package(pybind11 REQUIRED)
FILE(GLOB SRCFILES *.cpp)
pybind11_add_module(map_strings $(SRCFILES))

install (FILES map_strings.so DESTINATION lib)