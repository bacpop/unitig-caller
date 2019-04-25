/*
 * File: map_strings.cpp
 *
 * Reports exact matches of kmers to sequences
 *
 */

#include "map_strings.hpp"

void call_strings(std::vector<std::string>& assembly_list,
                  std::vector<std::string>& assembly_names,
                  std::vector<std::string>& query_list,
                  std::string& output_file,
                  size_t num_threads)
{
   assert(num_threads >= 1);

   // Read all sequences into memory as Fasta objects
   std::cerr << "Reading reference sequences into memory..." << std::endl;
   std::vector<Fasta> sequences;
   auto name_it = assembly_names.begin();
   for (auto file_it = assembly_list.begin(); file_it != assembly_list.end(); ++file_it)
   {
      sequences.push_back(Fasta(*name_it, *file_it));
      name_it++;
   }

   std::cerr << "Calling unitigs..." << std::endl;
   std::ofstream pres_ofs(output_file.c_str());

   // Create threaded queue for distance calculations
   std::queue<std::future<std::vector<std::string>>> map_threads;
   const unsigned long int calc_per_thread = (unsigned long int)sequences.size() / num_threads;
   const unsigned int num_big_threads = sequences.size() % num_threads;

   size_t start = 0;
   std::vector<size_t> start_points;
   for (unsigned int thread_idx = 0; thread_idx < threads; ++thread_idx) // Loop over threads
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

   // Run searches. Each thread looks at a chunk of the input fastas, then looping
   // over all unitig queries (and their reverse complements)
   for (auto unitig_it = query_list.begin(); unitig_it != query_list.end(); unitig_it++)
   {
      for (unsigned int thread_idx = 0; thread_idx < threads; ++thread_idx)
      {
         // Set the thread off
         map_threads.push(std::async(std::launch::async, seq_search,
                                                         std::cref(*unitig_it),
                                                         std::cref(sequences)
                                                         start_points[thread_idx],
                                                         start_points[thread_idx + 1]));
      }

      pres_ofs << *unitigs_it << " |";
      while (!map_threads.empty())
      {
         std::vector<std::string> present = distance_calculations.front().get();
         distance_calculations.pop();

         for (auto pres_it = present.begin() ; pres_it < present.end(); ++pres_it)
         {
            pres_ofs << " " << *pres_it << ":1";
         }
      }
      os << std::endl;
   }

   std::cerr << "Done." << std::endl;
}

std::vector<std:string> seq_search(std::string& query, std::vector<Fasta>& sequences, size_t start, size_t end)
{
   std::string rev_query = rev_comp(query);
   std::vector<std::string> present;
   for (auto fasta_it = sequences.begin() + start; fasta_it != sequences.begin() + end; fasta_it++)
   {
      if (fasta_it->hasSeq(query) || fasta_it->hasSeq(rev_query))
      {
         present.push_back(fasta_it->get_name());
      }
   }
   return present;
}

