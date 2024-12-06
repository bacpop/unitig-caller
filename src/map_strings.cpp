/*
 * File: map_strings.cpp
 *
 * Reports exact matches of kmers to sequences
 *
 */

// kseq headers are already included by bifrost
#include "map_strings.hpp"

// sdsl headers
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
typedef sdsl::csa_wt<> fm_index_coll;

// code from
// https://stackoverflow.com/questions/735204/convert-a-string-in-c-to-upper-case
char ascii_toupper_char(char c) {
  return ('a' <= c && c <= 'z')
             ? c ^ 0x20
             : c;
}

fm_index_coll index_fasta(const std::string &fasta_file,
                          const bool &write_idx) {
  fm_index_coll ref_index;
  // create fm index file name
  std::string idx_file_name = fasta_file + ".fm";
  // if fm_index not available, generate it
  if (!load_from_file(ref_index, idx_file_name)) {
    std::string reference_seq;
    // open the file handler
    gzFile fp = gzopen(fasta_file.c_str(), "r");
    if (fp == 0) {
      perror("fopen");
      exit(1);
    }
    // initialize seq
    kseq_t *seq = kseq_init(fp);
    // read sequence
    int l;
    while ((l = kseq_read(seq)) >= 0) {
      reference_seq += seq->seq.s;
      reference_seq += ",";
    }
    // destroy seq and fp objects
    kseq_destroy(seq);
    gzclose(fp);

    // Convert all to uppercase
    for (char &c : reference_seq) {
      c = ascii_toupper_char(c);
    }

    sdsl::construct_im(ref_index, reference_seq, 1); // generate index
    if (write_idx) {
      store_to_file(ref_index, idx_file_name); // save it
    }
  }
  return ref_index;
}

// search for a specific sequence within an fm index array
bool seq_search(std::string &query, const fm_index_coll &ref_idx) {
  // Convert query to uppercase
  for (char &c : query) {
    c = ascii_toupper_char(c);
  }

  bool present = false;
  // count number of occurrences in positive strand
  size_t query_count = sdsl::count(ref_idx, query.begin(), query.end());
  // if not found, check reverse strand
  if (query_count == 0) {
    // Revcomp from bifrost
    const std::string rev_query = reverse_complement(query);
    query_count = sdsl::count(ref_idx, rev_query.begin(), rev_query.end());
  }
  if (query_count != 0) {
    present = 1;
  }
  return present;
}

void call_strings(const std::vector<std::string> &assembly_list,
                  const std::vector<std::string> &assembly_names,
                  const std::vector<std::string> &query_list,
                  const std::string &output_file, const bool write_idx,
                  const size_t num_threads) {
  // Read all sequences into memory as Fasta objects (threaded)
  std::cerr << "Constructing indexes for all input sequences..." << std::endl;
  std::vector<fm_index_coll> seq_idx(assembly_list.size());
  #pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int file_idx = 0; file_idx < seq_idx.size(); ++file_idx) {
    seq_idx[file_idx] = index_fasta(assembly_list[file_idx], write_idx);
  }

  std::cerr << "Calling unitigs..." << std::endl;
  std::ofstream pres_ofs(output_file.c_str());

  // Run searches. Each thread looks at a chunk of the input fastas, then
  // looping over all unitig queries (and their reverse complements)
  #pragma omp parallel for schedule(dynamic) num_threads(num_threads)
  for (int unitig_idx = 0; unitig_idx < query_list.size(); ++unitig_idx) {
    std::vector<std::string> present;
    for (int fm_idx = 0; fm_idx < seq_idx.size(); ++fm_idx) {
      if (seq_search(query_list[unitig_idx], seq_idx[fm_idx])) {
        present.push_back(assembly_names[fm_idx]);
      }
    }
    // Print results if found
    if (present.size() > 0) {
      #pragma omp critical
      {
        pres_ofs << query_list[unitig_idx] << " |";
        for (auto pres_it = present.begin(); pres_it < present.end(); ++pres_it) {
          pres_ofs << " " << *pres_it << ":1";
        }
        pres_ofs << std::endl;
      }
    }
  }

  std::cerr << "Done." << std::endl;
}

