/*
 *
 * map_strings.hpp
 * Header file for map_strings
 *
 */

// C/C++/C++11/C++17 headers
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
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
#include <experimental/filesystem>

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// seqan3 headers
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/view/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <cereal/archives/binary.hpp>

namespace py = pybind11;
using namespace seqan3;

typedef fm_index<true, default_sdsl_index_type> fasta_fm_index;

// Constants
const std::string VERSION = "1.1.0";

// Structs

// Function headers
// map_strings.cpp
void call_strings(const std::vector<std::string>& assembly_list,
                  const std::vector<std::string>& assembly_names,
                  const std::vector<std::string>& query_list,
                  const std::string& output_file,
                  const bool write_idx = 1,
                  const size_t num_threads = 1);
std::vector<fasta_fm_index> index_fastas(const std::vector<std::string>& fasta_files,
                                        const size_t start,
                                        const size_t end,
                                        const bool write_idx = 1);
std::vector<std::string> seq_search(const dna5_vector& query,
                                    const std::vector<fasta_fm_index>& sequences,
                                    const std::vector<std::string>& names,
                                    const size_t start,
                                    const size_t end);

// map_bindings.cpp
int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    bool write_idx = 1,
                    size_t num_threads = 1);

