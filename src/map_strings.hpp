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
#include <queue>
#include <functional>
#include <future>
#include <thread>
#include <list>
#include <assert.h>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// seqan3 headers
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

namespace py = pybind11;
using namespace seqan3;

// Constants
const std::string VERSION = "0.3.0";

// Structs

// Function headers
// map_strings.cpp
void call_strings(const std::vector<std::string>& assembly_list,
                  const std::vector<std::string>& assembly_names,
                  const std::vector<std::string>& query_list,
                  const std::string& output_file,
                  const size_t num_threads = 1);
std::vector<std::string> seq_search(const std::string& query,
                                    const std::vector<fm_index>& sequences,
                                    const std::vector<std::string>& names,
                                    const size_t start,
                                    const size_t end);

// map_bindings.cpp
int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    size_t num_threads = 1);

