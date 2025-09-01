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

#ifdef FS_EXP
#include <experimental/filesystem>
#else
#include <filesystem>
#endif

// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Bifrost headers
#include <bifrost/ColoredCDBG.hpp>

namespace py = pybind11;

// return pair typedef
typedef std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>> ReturnPair;

// Function headers
// map_strings.cpp
void call_strings(const std::vector<std::string>& assembly_list,
                  const std::vector<std::string>& assembly_names,
                  std::vector<std::string>& query_list,
                  const std::string& output_file,
                  const bool write_idx = 1,
                  const size_t num_threads = 1);

// bifrost.cpp
// build graph from refs/reads
ColoredCDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool& is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool& write_graph,
                          const std::string& output_prefix);

// generate colours for a unitig
template <class T, class U, bool is_const>
std::vector<bool> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                   const size_t nb_colours,
                                   const size_t position,
                                   const bool single_mapping);

// negate colours from untigs
std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2);

// parse a fasta and return sequences
std::vector<std::string> parse_fasta (const std::string& fasta);

// query whether a unitig exists in a graph in it's entirity, if so return colour vector. If not, return empty vector.
void query_unitig (const ColoredCDBG<>& ccdbg, 
                    const std::vector<std::string>& query_list, 
                    const size_t& nb_colours, 
                    const std::string& out_path,
                    const std::vector<std::string>& input_colour_pref,
                    bool rtab,
                    bool pyseer);

// call unitigs and return their colours within a graph
void call_unitigs(const ColoredCDBG<>& ccdbg, 
                    const std::string& out_path,
                    const std::vector<std::string>& input_colour_pref,
                    bool rtab,
                    bool pyseer);

// map_bindings.cpp
int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    bool write_idx = 1,
                    size_t num_threads = 1);

std::vector<std::string> py_uc_exists(const std::string &graphfile,
                          const std::string &coloursfile, const bool call,
                          const std::string &query_file, const std::string &outpref,
                          size_t num_threads, bool rtab, bool pyseer) ;

std::vector<std::string> py_uc_build(const std::string &infile1, const int &kmer,
                       const bool call, const std::string &query_file,
                       size_t num_threads, bool is_ref, const bool write_graph, 
                       const std::string &outpref, bool rtab, bool pyseer,
                       const std::string &infile2) ;
