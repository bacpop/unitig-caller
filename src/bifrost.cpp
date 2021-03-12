#include "map_strings.hpp"

// build graph from refs/reads
ColoredCDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool& is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool& write_graph,
                          const std::string& output_prefix)
{
    std::ifstream infile1(infile_1);
    std::ifstream infile2(infile_2);
    CCDBG_Build_opt opt;


    opt.k = kmer;
    opt.nb_threads = threads;
    opt.verbose = verb;
    opt.prefixFilenameOut = output_prefix;

    std::string filename;
    if (is_ref && (infile_2 == "NA")) {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
    } else if (!is_ref && (infile_2 == "NA"))
    {
        while (std::getline(infile1, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    } else {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
        while (std::getline(infile2, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    }

    ColoredCDBG<> ccdbg(opt.k);
    ccdbg.buildGraph(opt);
    ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
    ccdbg.buildColors(opt);

    if (write_graph)
    {
        ccdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
    }

    return ccdbg;
}

// generate colours for a unitig
template <class T, class U, bool is_const>
std::vector<bool> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                   const size_t nb_colours,
                                   const size_t position,
                                   const bool single_mapping)
{
    // get colours information for unitig
    const auto colourset = unitig.getData()->getUnitigColors(unitig);
    std::vector<bool> colours_arr(nb_colours, 0);

    // initialise a iterator, will only determine colours of single kmer, have to do for head and tail as may be different
    UnitigColors::const_iterator it_uc = colourset->begin(unitig);
    UnitigColors::const_iterator it_uc_end = colourset->end();

    // if not single mapping, look at specific position
    if (!single_mapping)
    {
        // iterate over specified kmer colours, update colours array with presence of those colours
        for (it_uc; it_uc != it_uc_end; it_uc++)
        {
            if (it_uc.getKmerPosition() == position)
            {
                colours_arr[it_uc.getColorID()] = 1;
            }
        }
    }
    // else, look at only position in unitig given by mapping, ignore position
    else {
        // iterate over specified kmer colours, update colours array with presence of those colours
        for (it_uc; it_uc != it_uc_end; it_uc++)
        {
            colours_arr[it_uc.getColorID()] = 1;
        }
    }
    return colours_arr;
}

// negate colour arrays
std::vector<bool> negate_colours_array(const std::vector<bool>& array1, const std::vector<bool>& array2)
{
    std::vector<bool> output_array = array1;
    for (size_t i = 0; i < array1.size(); i++)
    {
        if (array1[i] == 1 && array2[i] == 0)
        {
            output_array[i] = 0;
        }
    }
    return output_array;
}

// parse a fasta and return sequences
std::vector<std::string> parse_fasta (const std::string& fasta)
{
    std::vector<std::string> seq_list;
    std::string line, DNA_sequence;
    std::ifstream infile(fasta);
    // parse fasta
    while (std::getline(infile, line)) {
        // remove new line characters
        line.erase(std::remove(line.begin(), line.end(), '\n'),
                    line.end());

        // erase DNA_sequence
        DNA_sequence.clear();

        // line may be empty so you *must* ignore blank lines
        if (line.empty())
            continue;

        if (line[0] != '>') {
            DNA_sequence = line;
        }

        // add to seq list
        if (DNA_sequence.size() != 0)
            seq_list.push_back(DNA_sequence);
    }

    return seq_list;
}

// query whether a unitig exists in a graph in it's entirity, if so return colour vector. If not, return empty vector.
std::vector<bool> query_unitig (const ColoredCDBG<>& ccdbg, const std::string& query, const size_t& nb_colours)
{
    //split query into sequence of kmers
    const char *query_str = query.c_str();

    // find first kmer within query in graph
    KmerIterator it_km(query_str), it_km_end;

    // query the kmer in the graph, find the first unitig
    auto unitig_current = ccdbg.find(it_km->first);

    // make a const copy of the first unitig map to enable colour generation
    const auto unitig_first = unitig_current;

    // get the colours for the unitig
    std::vector<bool> query_colours;
    std::vector<bool> query_colours_head;
    std::vector<bool> query_colours_tail;

//    // initialise check if kmer is first or last in sequence
//    bool is_first = true;
//    bool is_last = false;

    // initialise unitig_found check
    bool unitig_present = true;

    // iterate to next kmer iterator
    it_km++;

    // iterate over forward strand of unitig
    if (!unitig_current.isEmpty)
    {
        for (it_km; it_km != it_km_end; it_km++)
        {
            // query the kmer in the graph
            auto unitig_current = ccdbg.find(it_km->first);

            // check if unitig_current is true, if not then break and set unitig_present to false
            if (unitig_current.isEmpty)
            {
                unitig_present = false;
                break;
            }

            // else, continue, the next iteration unitig_current will be set
        }
    } else {
        unitig_present = false;
    }

    // if unitig not present, clear the colours
    if (unitig_present)
    {
        // get head and tail colours
        query_colours_head = generate_colours(unitig_first, nb_colours, 0, true);
        query_colours_tail = generate_colours(unitig_current, nb_colours, 0, true);

        // if colours are not the same, negate head and tail colours
        if (query_colours_head != query_colours_tail) {
            query_colours = std::move(negate_colours_array(query_colours_head, query_colours_tail));
        } else{
            query_colours = std::move(query_colours_head);
        }
    }

    return query_colours;
}

// call unitigs and return their colours within a graph
std::unordered_map<std::string, std::vector<bool>> call_unitigs(const ColoredCDBG<>& ccdbg)
{
    // initialise unitig map to return
    std::unordered_map<std::string, std::vector<bool>> unitig_map;

//    // get kmer size
//    const int kmer = ccdbg.getK();

    // get the number of colours
    const size_t nb_colours = ccdbg.getNbColors();

    for (const auto um : ccdbg)
    {
        // get tail kmer position as len of unitig in kmers minus 1 (zero basedness)
        const int tail_pos = um.len - 1;

        // generate colours for unitig
        std::vector<bool> um_colours_head = generate_colours(um, nb_colours, 0, false);
        std::vector<bool> um_colours_tail = generate_colours(um, nb_colours, tail_pos, false);

        // negate colours if head and tail aren't equal, otherwise set to head
        std::vector<bool> um_colours;
        if (um_colours_head != um_colours_tail) {
            um_colours = std::move(negate_colours_array(um_colours_head, um_colours_tail));
        } else{
            um_colours = std::move (um_colours_head);
        }

        // generate string for unitig
        std::string um_seq = um.referenceUnitigToString();

        // append to unitig_map
        unitig_map.insert(std::make_pair(um_seq, um_colours));
    }

    return unitig_map;
}


void print_unitigs (const std::pair<std::unordered_map<std::string, std::vector<bool>>, std::vector<std::string>>& return_pair,
        const std::string& outfile_name)
{
    ofstream outfile;
    outfile.open(outfile_name);

    cout << "Printing unitigs..." << endl;

    for (const auto& unitig : return_pair.first)
    {
        // generate string for colours
        std::string colours;
        for (const auto& i : unitig.second)
        {
            colours += std::to_string(i);
        }
        // append to file
        outfile << unitig.first << "\t" << colours << "\n";
    }

    outfile.close();
}