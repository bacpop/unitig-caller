/*
 * File: fasta.cpp
 *
 * Helper functions for the fasta sequence class
 *
 */

#include "fasta.hpp"

// Read in sequences to memory from file
Fasta::Fasta(const std::string& obj_name, const std::string& filename)
   :name(obj_name)
{
   std::ifstream ist(filename.c_str());

   if (!ist)
   {
      throw std::runtime_error("Could not open fasta file " + filename + "\n");
   }

   std::string line_in;
   std::string seq, contig_name = "";

   std::getline(ist, line_in);

   if (ist.eof())
   {
       throw std::runtime_error("Could not read fasta file " + filename + "\n");
   }
   // Read header line, which starts with >
   else if (line_in[0] == '>')
   {
      contig_name = line_in.substr(1);
   }
   else
   {
      throw std::runtime_error("Error reading fasta file " + filename + "\n"
         "Header should start with '>', but has:\n" + line_in + "\n");
   }

   while (!ist.eof())
   {
      // New contig, new sequence
      if (ist.peek() == '>')
      {
         sequence_names.push_back(contig_name);
         sequences.push_back(seq);

         seq = "";
         std::getline(ist, line_in);
         contig_name = line_in.substr(1);
      }
      else
      {
         // Concatenate sequence lines
         std::getline(ist, line_in);
         seq += line_in;
      }
   }

   // Add in final contig
   sequence_names.push_back(contig_name);
   sequences.push_back(seq);

}

// Simple presence/absence test
bool Fasta::hasSeq(const std::string& search)
{
   bool found = 0;

   // Search each sequence in fasta
   for (std::vector<std::string>::iterator it = Fasta::sequences.begin(); it != Fasta::sequences.end(); ++it)
   {
      size_t hit_position = ;
      if (it->find(search) != std::string::npos)
      {
         found = 1;
         break;
      }
   }

   return found;
}

std::string rev_comp(std::string& seq)
{
   std::string rev_seq(seq);
   auto lambda = [](const char c)
   {
      switch (c)
      {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            throw std::runtime_error("Invalid nucleotide");
        }
    };

    std::transform(rev_seq.cbegin(), rev_seq.cend(), rev_seq.begin(), lambda);
    std::reverse(rev_seq.begin(), rev_seq.end());
    return rev_seq;
}