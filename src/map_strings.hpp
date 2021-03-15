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
#include <seqan3/search/fm_index/all.hpp>
#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/search/search.hpp>
//#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

// Bifrost headers
#include <bifrost/ColoredCDBG.hpp>

namespace py = pybind11;
//using namespace seqan3;

// fmindex typedef
using cust_sdsl_wt_index_type = sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
        sdsl::rank_support_v<>,
        sdsl::select_support_scan<>,
        sdsl::select_support_scan<0>>,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::plain_byte_alphabet>;
typedef seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection, cust_sdsl_wt_index_type> fasta_fm_index;
using seqan3::operator""_dna5;

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
std::vector<std::string> seq_search(const seqan3::dna5_vector& query,
                                    const std::vector<fasta_fm_index>& seq_idx,
                                    const std::vector<std::string>& names,
                                    const size_t start,
                                    const size_t end);

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
std::vector<bool> query_unitig (const ColoredCDBG<>& ccdbg, const std::string& query, const size_t& nb_colours);

// call unitigs and return their colours within a graph
std::unordered_map<std::string, std::vector<bool>> call_unitigs(const ColoredCDBG<>& ccdbg);

void print_unitigs (const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>>& return_pair,
                    const std::string& outfile_name);

// map_bindings.cpp
int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    bool write_idx = 1,
                    size_t num_threads = 1);

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_call_exists (const std::string& graphfile,
                                                     const std::string& coloursfile,
                                                     size_t num_threads);

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_call_build (const std::string& infile1,
                                                    const int& kmer,
                                                    size_t num_threads,
                                                    bool is_ref,
                                                    const bool write_graph,
                                                    const std::string& infile2);

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_query_exists (const std::string& graphfile,
                                                      const std::string& coloursfile,
                                                      const std::string& query_file,
                                                      size_t num_threads);

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_query_build (const std::string& infile1,
                                                     const int& kmer,
                                                     const std::string& query_file,
                                                     size_t num_threads,
                                                     bool is_ref,
                                                     const bool write_graph,
                                                     const std::string& infile2);
