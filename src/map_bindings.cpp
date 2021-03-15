/*
 * map_bindings.cpp
 * Python bindings for map_strings
 *
 */

#include "map_strings.hpp"

int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file,
                    bool write_idx,
                    size_t num_threads)
{
    // Check input

    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // Convert python objs to C++
    // Here done automatically with pybind11/stl.h

    // call pure C++ function
    call_strings(assembly_list, assembly_names, query_list, output_file, write_idx, num_threads);

    // return success
    return 1;
}

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_call_exists (const std::string& graphfile,
                                                     const std::string& coloursfile,
                                                     size_t num_threads) {
    // Set number of threads
    if (num_threads < 1) {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Reading coloured compacted DBG..." << endl;

    // read in graph
    ColoredCDBG<> ccdbg;
    ccdbg.read(graphfile, coloursfile, num_threads);

    cout << "Calling unitigs within population..." << endl;

    // get colour names
    std::vector<std::string> input_colour_files = ccdbg.getColorNames();
    std::vector<std::string> input_colour_pref;

    // generate file prefixes for colours
    for (const auto& file : input_colour_files)
    {
        std::filesystem::path p(file);
        input_colour_pref.push_back(p.stem());
    }

    auto unitig_map = call_unitigs(ccdbg);
    const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>> return_pair = std::make_pair(unitig_map, input_colour_pref);

    return return_pair;
}

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_call_build (const std::string& infile1,
                                                    const int& kmer,
                                                    size_t num_threads,
                                                    bool is_ref,
                                                    const bool write_graph,
                                                    const std::string& infile2) {
    // Set number of threads
    if (num_threads < 1) {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Building coloured compacted DBG..." << endl;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    ColoredCDBG<> ccdbg;

    // build graph, write graph if specified
    size_t lastindex = infile1.find_last_of(".");
    std::string outgraph = infile1.substr(0, lastindex);
    ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

    cout << "Calling unitigs within population..." << endl;

    // get colour names
    std::vector<std::string> input_colour_files = ccdbg.getColorNames();
    std::vector<std::string> input_colour_pref;

    // generate file prefixes for colours
    for (const auto& file : input_colour_files)
    {
        std::filesystem::path p(file);
        input_colour_pref.push_back(p.stem());
    }

    auto unitig_map = call_unitigs(ccdbg);
    const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>> return_pair = std::make_pair(unitig_map, input_colour_pref);

    return return_pair;
}

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_query_exists (const std::string& graphfile,
                                                      const std::string& coloursfile,
                                                      const std::string& query_file,
                                                      size_t num_threads) {
    // Set number of threads
    if (num_threads < 1) {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Reading coloured compacted DBG..." << endl;

    // read in graph
    ColoredCDBG<> ccdbg;
    ccdbg.read(graphfile, coloursfile, num_threads);

    std::vector<std::string> query_list = parse_fasta(query_file);

    // get kmer size
    const int kmer = ccdbg.getK();

    // get the number of colours
    const size_t nb_colours = ccdbg.getNbColors();

    std::unordered_map<std::string, std::vector<bool>> query_colours_map;

    cout << "Querying unitigs within population..." << endl;
    for (const auto& query : query_list)
    {
        // run query of colours
        std::vector<bool> query_colours = query_unitig(ccdbg, query, nb_colours);

        // Add colours to map. If unitig not found, query colours will be empty
        query_colours_map[query] = std::move(query_colours);
    }

    // get colour names
    std::vector<std::string> input_colour_files = ccdbg.getColorNames();
    std::vector<std::string> input_colour_pref;

    // generate file prefixes for colours
    for (const auto& file : input_colour_files)
    {
        std::filesystem::path p(file);
        input_colour_pref.push_back(p.stem());
    }

    const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>> return_pair = std::make_pair(query_colours_map, input_colour_pref);

    return return_pair;
}

std::pair<std::unordered_map<std::string, std::vector<bool>>,
        std::vector<std::string>> py_uc_query_build (const std::string& infile1,
                                                     const int& kmer,
                                                     const std::string& query_file,
                                                     size_t num_threads,
                                                     bool is_ref,
                                                     const bool write_graph,
                                                     const std::string& infile2) {
    // Set number of threads
    if (num_threads < 1) {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Building coloured compacted DBG..." << endl;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    ColoredCDBG<> ccdbg;

    // build graph, write graph if specified
    size_t lastindex = infile1.find_last_of(".");
    std::string outgraph = infile1.substr(0, lastindex);
    ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

    std::vector<std::string> query_list = parse_fasta(query_file);

    // get the number of colours
    const size_t nb_colours = ccdbg.getNbColors();

    std::unordered_map<std::string, std::vector<bool>> query_colours_map;

    cout << "Querying unitigs within population..." << endl;
    for (const auto& query : query_list)
    {
        // run query of unitigs
        std::vector<bool> query_colours = query_unitig(ccdbg, query, nb_colours);

        // Add colours to map. If unitig not found, query colours will be empty
        query_colours_map[query] = std::move(query_colours);
    }

    // get colour names
    std::vector<std::string> input_colour_files = ccdbg.getColorNames();
    std::vector<std::string> input_colour_pref;

    // generate file prefixes for colours
    for (const auto& file : input_colour_files)
    {
        std::filesystem::path p(file);
        input_colour_pref.push_back(p.stem());
    }

    const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>> return_pair = std::make_pair(query_colours_map, input_colour_pref);

    return return_pair;
}

PYBIND11_MODULE(map_strings, m)
{
  m.doc() = "Finds presence/absence of substrings";

  m.def("call", &py_call_strings, "Print presence absence to file",
        py::arg("assembly_files"),
        py::arg("assembly_names"),
        py::arg("queries"),
        py::arg("output_file"),
        py::arg("write_idx") = 1,
        py::arg("threads") = 1);

  m.def("call_unitigs_existing", &py_uc_call_exists, "Calls unitigs and their colours in an existing Bifrost graph",
          py::arg("graphfile"),
          py::arg("coloursfile"),
          py::arg("num_threads") = 1);

  m.def("call_unitigs_build", &py_uc_call_build, "Builds and then calls unitigs in Bifrost graph",
          py::arg("infile1"),
          py::arg("kmer"),
          py::arg("num_threads") = 1,
          py::arg("is_ref") = 1,
          py::arg("write_graph") = 0,
          py::arg("infile2") = "NA");

  m.def("query_unitigs_existing", &py_uc_query_exists, "Queries unitigs and their colours in an existing Bifrost graph",
          py::arg("graphfile"),
          py::arg("coloursfile"),
          py::arg("query_file"),
          py::arg("num_threads") = 1);

  m.def("query_unitigs_build", &py_uc_query_build, "Builds and then queries unitigs in Bifrost graph",
          py::arg("infile1"),
          py::arg("kmer"),
          py::arg("query_file"),
          py::arg("num_threads") = 1,
          py::arg("is_ref") = 1,
          py::arg("write_graph") = 0,
          py::arg("infile2") = "NA");
}
