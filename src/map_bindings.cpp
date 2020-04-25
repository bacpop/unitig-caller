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
}