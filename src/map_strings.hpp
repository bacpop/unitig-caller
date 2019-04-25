/*
 *
 * map_strings.hpp
 * Header file for map_strings
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <functional>
#include <future>
#include <thread>
#include <list>
#include <assert.h>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Classes
#include "fasta.hpp"

// Constants
const std::string VERSION = "0.1.0";

// Structs

// Function headers
// map_strings.cpp
void call_strings(std::vector<std::string>& assembly_list,
                  std::vector<std::string>& assembly_names,
                  std::vector<std::string>& query_list,
                  std::string& output_file,
                  size_t num_threads = 1);
std::vector<std:string> seq_search(std::string& query, std::vector<Fasta>& sequences, size_t start, size_t end);

// map_bindings.cpp
int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    size_t num_threads = 1);

