/*
 * fasta.hpp
 * Header file for fasta class
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <mutex>

// Contains sequences and contig names
class Fasta
{
   public:
      // Initialisation
      Fasta(const std::string& obj_name, const std::string& file);

      // Complex operations
      bool hasSeq(const std::string& search);

      // nonmodifying operations
      std::string get_name() const { return name; }
      static bool compareFasta(Fasta lhs, Fasta rhs) { return (lhs.name < rhs.name); }

   private:
      // Sequences and sequence names necessarily in same order as read in
      std::vector<std::string> sequences;
      std::vector<std::string> sequence_names;

      std::string name;
};

std::string rev_comp(std::string& seq);