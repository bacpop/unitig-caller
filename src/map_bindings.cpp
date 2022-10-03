/*
 * map_bindings.cpp
 * Python bindings for map_strings
 *
 */

#include "map_strings.hpp"

int py_call_strings(std::vector<std::string> assembly_list,
                    std::vector<std::string> assembly_names,
                    std::vector<std::string> query_list,
                    std::string output_file, bool write_idx,
                    size_t num_threads) {
  // Check input

  // Set number of threads
  if (num_threads < 1) {
    num_threads = 1;
  }

  // Convert python objs to C++
  // Here done automatically with pybind11/stl.h

  // call pure C++ function
  call_strings(assembly_list, assembly_names, query_list, output_file,
               write_idx, num_threads);

  // return success
  return 1;
}

ReturnPair py_uc_exists(const std::string &graphfile,
                        const std::string &coloursfile, const bool call,
                        const std::string &query_file, size_t num_threads) {
  // Set number of threads
  if (num_threads < 1) {
    num_threads = 1;
  }

  // read in compact coloured DBG
  cout << "Reading coloured compacted DBG..." << endl;

  // read in graph
  ColoredCDBG<> ccdbg;
  ccdbg.read(graphfile, coloursfile, num_threads);

  // get colour names
  std::vector<std::string> input_colour_files = ccdbg.getColorNames();
  std::vector<std::string> input_colour_pref;

  // generate file prefixes for colours
  for (const auto &file : input_colour_files) {
#ifdef FS_EXP
    std::experimental::filesystem::path p(file);
#else
    std::filesystem::path p(file);
#endif
    input_colour_pref.push_back(p.stem());
  }

  std::unordered_map<std::string, std::vector<bool>> unitig_map;

  // get the number of colours
  const size_t nb_colours = ccdbg.getNbColors();

  if (call) {
    cout << "Calling unitigs within population..." << endl;
    unitig_map = call_unitigs(ccdbg);
  } else {
    cout << "Querying unitigs within population..." << endl;

    std::vector<std::string> query_list = parse_fasta(query_file);

    for (const auto &query : query_list) {
      // run query of colours
      std::vector<bool> query_colours = query_unitig(ccdbg, query, nb_colours);

      // Add colours to map. If unitig not found, query colours will be empty
      unitig_map[query] = std::move(query_colours);
    }
  }

  const ReturnPair return_pair = std::make_pair(unitig_map, input_colour_pref);

  return return_pair;
}

ReturnPair py_uc_build(const std::string &infile1, const int &kmer,
                       const bool call, const std::string &query_file,
                       size_t num_threads, bool is_ref, const bool write_graph,
                       const std::string &infile2) {
  // Set number of threads
  if (num_threads < 1) {
    num_threads = 1;
  }

  // read in compact coloured DBG
  cout << "Building coloured compacted DBG..." << endl;

  if (infile2 != "NA") {
    is_ref = 0;
  }

  // build graph, write graph if specified
  size_t lastindex = infile1.find_last_of(".");
  std::string outgraph = infile1.substr(0, lastindex);
  ColoredCDBG<> ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads,
                                   false, write_graph, outgraph);

  // get colour names
  std::vector<std::string> input_colour_files = ccdbg.getColorNames();
  std::vector<std::string> input_colour_pref;

  // generate file prefixes for colours
  for (const auto &file : input_colour_files) {
#ifdef FS_EXP
    std::experimental::filesystem::path p(file);
#else
    std::filesystem::path p(file);
#endif
    input_colour_pref.push_back(p.stem());
  }

  std::unordered_map<std::string, std::vector<bool>> unitig_map;

  // get the number of colours
  const size_t nb_colours = ccdbg.getNbColors();

  if (call) {
    cout << "Calling unitigs within population..." << endl;
    unitig_map = call_unitigs(ccdbg);
  } else {
    cout << "Querying unitigs within population..." << endl;

    std::vector<std::string> query_list = parse_fasta(query_file);

    for (const auto &query : query_list) {
      // run query of colours
      std::vector<bool> query_colours = query_unitig(ccdbg, query, nb_colours);

      // Add colours to map. If unitig not found, query colours will be empty
      unitig_map[query] = std::move(query_colours);
    }
  }

  const ReturnPair return_pair = std::make_pair(unitig_map, input_colour_pref);

  return return_pair;
}

PYBIND11_MODULE(unitig_query, m) {
  m.doc() = "Finds presence/absence of substrings";

  m.def("call", &py_call_strings, "Print presence absence to file",
        py::arg("assembly_files"), py::arg("assembly_names"),
        py::arg("queries"), py::arg("output_file"), py::arg("write_idx") = 1,
        py::arg("threads") = 1);

  m.def("call_unitigs_existing", &py_uc_exists,
        "Call/queries unitigs and their colours in an existing Bifrost graph",
        py::arg("graphfile"), py::arg("coloursfile"), py::arg("call"),
        py::arg("query_file"), py::arg("num_threads") = 1);

  m.def("call_unitigs_build", &py_uc_build,
        "Builds and then calls/queries unitigs in Bifrost graph",
        py::arg("infile1"), py::arg("kmer"), py::arg("call"),
        py::arg("query_file"), py::arg("num_threads") = 1,
        py::arg("is_ref") = 1, py::arg("write_graph") = 0,
        py::arg("infile2") = "NA");

  m.attr("version") = VERSION_INFO;
}
