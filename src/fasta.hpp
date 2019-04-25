/*
 * fasta.hpp
 * Header file for fasta class
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>

// Contains sequences and contig names
class Fasta
{
   public:
      // Initialisation
      Fasta(const std::string& obj_name, const std::string& file);

      // nonmodifying operations
      bool hasSeq(const std::string& search) const;
      std::string get_name() const { return name; }
      static bool compareFasta(Fasta lhs, Fasta rhs) { return (lhs.name < rhs.name); }

   private:
      // Sequences and sequence names necessarily in same order as read in
      std::vector<std::string> sequences;
      std::vector<std::string> sequence_names;

      std::string name;
};

std::string rev_comp(const std::string& seq);