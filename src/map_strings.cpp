/*
 * File: map_strings.cpp
 *
 * Reports exact matches of kmers to sequences
 *
 */

#include "map_strings.hpp"

void call_strings(const std::vector<std::string>& assembly_list,
                  const std::vector<std::string>& assembly_names,
                  const std::vector<std::string>& query_list,
                  const std::string& output_file,
                  const bool write_idx,
                  const size_t num_threads)
{
   // Create threaded queue for computation
   assert(num_threads >= 1);
   const unsigned long int calc_per_thread = (unsigned long int)assembly_list.size() / num_threads;
   const unsigned int num_big_threads = assembly_list.size() % num_threads;

   size_t start = 0;
   std::vector<size_t> start_points;
   for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) // Loop over threads
   {
      start_points.push_back(start);

      // First 'big' threads have an extra job
      unsigned long int thread_jobs;
      if (thread_idx < num_big_threads)
      {
         start += calc_per_thread + 1;
      }
      else
      {
         start += calc_per_thread;
      }
   }
   start_points.push_back(start);

   // Read all sequences into memory as Fasta objects (threaded)
   std::cerr << "Constructing indexes for all input sequences..." << std::endl;
   std::queue<std::future<std::vector<fasta_fm_index>>> index_threads;
   for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx)
   {
      // Set the thread off
      index_threads.push(std::async(std::launch::async, index_fastas,
                                                        std::cref(assembly_list),
                                                        start_points[thread_idx],
                                                        start_points[thread_idx + 1],
                                                        write_idx));
   }

   // Get results from thread
   std::vector<fasta_fm_index> seq_idx;
   while (!index_threads.empty())
   {
      std::vector<fasta_fm_index> returned_indexes = index_threads.front().get();
      index_threads.pop();

      seq_idx.insert(seq_idx.end(), returned_indexes.begin(), returned_indexes.end());
   }

   std::cerr << "Calling unitigs..." << std::endl;
   std::ofstream pres_ofs(output_file.c_str());

   // Run searches. Each thread looks at a chunk of the input fastas, then looping
   // over all unitig queries (and their reverse complements)
   std::queue<std::future<std::vector<std::string>>> map_threads;
   for (auto unitig_it = query_list.begin(); unitig_it != query_list.end(); unitig_it++)
   {
      // debug_stream << *unitig_it << std::endl;
      seqan3::dna5_vector query = *unitig_it | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;

      for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx)
      {
         // Set the thread off
         map_threads.push(std::async(std::launch::async, seq_search,
                                                         std::cref(query),
                                                         std::cref(seq_idx),
                                                         std::cref(assembly_names),
                                                         start_points[thread_idx],
                                                         start_points[thread_idx + 1]));
      }

      // Get results from thread
      std::vector<std::string> present;
      while (!map_threads.empty())
      {
         std::vector<std::string> thread_present = map_threads.front().get();
         map_threads.pop();

         present.insert(present.end(), thread_present.begin(), thread_present.end());
      }

      // Print results if found
      if (present.size() > 0)
      {
         pres_ofs << *unitig_it << " |";
         for (auto pres_it = present.begin() ; pres_it < present.end(); ++pres_it)
         {
            pres_ofs << " " << *pres_it << ":1";
         }
         pres_ofs << std::endl;
      }
   }

   std::cerr << "Done." << std::endl;
}

std::vector<fasta_fm_index> index_fastas(const std::vector<std::string>& fasta_files,
                                        const size_t start,
                                        const size_t end,
                                        const bool write_idx)
{
   std::vector<fasta_fm_index> seq_idx;
   for (auto file_it = fasta_files.begin() + start; file_it != fasta_files.begin() + end; ++file_it)
   {
      std::string idx_file_name = *file_it + ".fm";

      // Read index if it already exists
      fasta_fm_index ref_index;
      if (std::filesystem::exists(idx_file_name))
      {
         {
            std::ifstream is{idx_file_name, std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(ref_index);
         }
      }
      else
      {
         // Create index
         seqan3::sequence_file_input reference_in{*file_it};
         std::vector<seqan3::dna5_vector> reference_seq;
         for (auto & [seq, id, qual] : reference_in)
         {
            reference_seq.push_back(std::move(seq));
         }
         ref_index = fasta_fm_index{reference_seq};

         // Write index to file
         if (write_idx)
         {
            std::ofstream os{idx_file_name, std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(ref_index);
         }
      }
      seq_idx.push_back(ref_index);
   }
   return seq_idx;
}

std::vector<std::string> seq_search(const seqan3::dna5_vector& query,
                                    const std::vector<fasta_fm_index>& seq_idx,
                                    const std::vector<std::string>& names,
                                    const size_t start,
                                    const size_t end)
{
   std::vector<std::string> present;
   auto name_it = names.begin() + start;
   for (auto ref_it = seq_idx.begin() + start; ref_it != seq_idx.begin() + end; ref_it++)
   {
      // debug_stream << *name_it << std::endl;

      int found = 0;
      auto results = search(query, *ref_it);
      found = (int)std::ranges::distance(results);
      // debug_stream << "There are " << results.size() << " hits.\n";
      // debug_stream << results << '\n';
      if (!found)
      {
         auto results = search(query | std::views::reverse | seqan3::views::complement, *ref_it);
         found = (int)std::ranges::distance(results);
         // debug_stream << "There are " << results.size() << " hits.\n";
         // debug_stream << results << '\n';
      }

      if (found)
      {
         // debug_stream << "found" << std::endl;
         present.push_back(*name_it);
      }
      name_it++;
   }
   return present;
}
