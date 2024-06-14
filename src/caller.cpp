#include <algorithm>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include "spoa/spoa.hpp"
#include "spoa/alignment_engine.hpp"
using namespace std;

#define PACKAGE_VERSION "0.1.6"

typedef struct compare_string_size_lt_t {
    const bool operator()(const std::string& first, const std::string& second) {
        return first.size() < second.size() || (first.size() == second.size() && first < second);
    }
} compare_string_size_lt_t;

typedef struct compare_string_size_gt_t {
    bool operator()(const std::string& first, const std::string& second) {
        return second.size() < first.size() || (second.size() == first.size() && second < first);
    }
} compare_string_size_gt_t;

static struct option options[] = {
    {"input", required_argument, 0, 'i'},
    {"match", required_argument, 0, 'A'},
    {"mismatch", required_argument, 0, 'B'},
    {"gap", required_argument, 0, 'O'},
    {"algorithm", required_argument, 0, 'a'},
    {"resort", no_argument, 0, 'r'},
    {"msa", no_argument, 0, 'm'},
    {"coverage", no_argument, 0, 'c'},
    {"left-align", no_argument, 0, 'l'},
    {"pairwise-msa", no_argument, 0, 'p'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

typedef struct consensus_opt_t {
    std::string input;
    int8_t match;
    int8_t mismatch;
    int8_t gap;
    int8_t algorithm;
    int8_t resort;
    bool coverage;
    bool msa;
    bool msa_extra;
    bool left_align;
    bool pairwise_msa;
    // for resort
    compare_string_size_lt_t _compare_lt;
    compare_string_size_gt_t _compare_gt;
} consensus_opt_t;

consensus_opt_t *consensus_opt_init()
{
    consensus_opt_t *opt = (consensus_opt_t*)calloc(1, sizeof(consensus_opt_t));
    opt->match           = 5;
    opt->mismatch        = -4;
    opt->gap             = -8;
    opt->algorithm       = 0;
    opt->resort          = 0;
    opt->coverage        = false;
    opt->msa             = false;
    opt->msa_extra       = false;
    opt->left_align      = false;
    opt->pairwise_msa    = false;
    return opt;
}

void left_align(std::string &consensus, std::string &sequence) {
    int msa_size = consensus.size();

    int left = 0;
    while(left < msa_size) {
        // find the left-moset '-' in a run of '-'
        if (sequence[left] != '-') {
            left++;
            continue;
        }
        int right = left;
        // move to base after the deletion
        while (right < msa_size && sequence[right] == '-') {
            right++;
        }
        if (msa_size == right) { // no more bases
            left = right;
            continue;
        }
        if (left == right) {
            left++;
            continue;
        }
        // examine base-by-base
        while (right < msa_size && consensus[left] == sequence[right]) {
            sequence[left] = sequence[right];
            sequence[right] = '-';
            left++;
            right++;
        }
        left = right;
    }
}

bool msa_all_dashes(const std::vector<std::string> &msa, const int index) {
    int i;
    for (i = 0; i < msa.size() && msa[i][index] == '-'; i++) {}  
    return (i == msa.size());
}

void left_align_msa(std::vector<std::string> &msa) {
    std::string consensus_msa = msa[msa.size()-1];

    // left-align
    for (int i = 0; i < msa.size()-1; i++) {
        left_align(consensus_msa, msa[i]);
    }

    // remove any positions in the MSA that are all '-'
    int left = 0;
    while (left < msa[0].size()) {
        // check if all positions are '-'
        int right = left;
        while (right < msa[0].size() && msa_all_dashes(msa, right)) {
            right++;
        }
        if (left == right) {
            left++;
            continue;
        }
        right--;
        int next_left = right + 1;

        // shift down
        while (right < msa[0].size()-1) {
            for (int j = 0; j < msa.size(); j++) {
                msa[j][left]  = msa[j][right+1];
                msa[j][right+1] = '-';
            }
            left++;
            right++;
        }

        // resize
        for (int j = 0; j < msa.size(); j++) {
            msa[j].resize(msa[j].size() - (right-left+1));
        }

        left = next_left;
    }
}

void process(std::unique_ptr<spoa::AlignmentEngine> &alignment_engine, std::string &name, std::vector<std::string> &sequences, consensus_opt_t *opt) {
    alignment_engine->Prealloc(sequences.size(), 4);


    switch(opt->resort) {
        case 0: break; // do nothing
        case 1: 
                std::sort(sequences.begin(), sequences.end(), opt->_compare_lt); 
                break;
        case 2: 
                compare_string_size_gt_t compare_gt;
                std::sort(sequences.begin(), sequences.end(), opt->_compare_gt); 
                break;
        default: 
                fprintf(stderr, "Bug: resort value '%d' not valid!\n", opt->resort);
                exit(1);
    }

    // create the graph for alignment
    //auto graph = spoa::createGraph();
    spoa::Graph graph{};

    // add the alignments to the graph
    for (const auto& it: sequences) {
        std::int32_t score = 0;
        auto alignment = alignment_engine->Align(it, graph, &score);
        graph.AddAlignment(alignment, it);
    }

    std::string consensus;
    if (opt->pairwise_msa) {
        // generate the consenuss
        consensus = graph.GenerateConsensus();
        // add the consenuss to the list of sequences
        sequences.insert(sequences.begin(), consensus);
        // create a enw graph for the second round of alignment
        graph = spoa::Graph();
        // add the alignments to the graph, including the consensus!
        for (const auto& it: sequences) {
            auto alignment = alignment_engine->Align(it, graph);
            graph.AddAlignment(alignment, it);
        }
    }

    // call the consensus
    if (opt->coverage) {
        std::vector<uint32_t> coverage;
        std::string consensus = graph.GenerateConsensus(&coverage, true);
        // reduce coverage by one if we added the consensus to the graph
        if (opt->pairwise_msa) {
            for (int i = 0; i < consensus.size(); i++) {
                uint8_t code = graph.coder(consensus[i]);
                coverage[code * consensus.size() + i]--;
            }
        }
        // coverage for each possible base
        fprintf(stdout, "%s\n%s\n", name.c_str(), consensus.c_str());
        for (uint32_t i = 0; i < graph.num_codes(); ++i) {
            fputc(graph.decoder(i), stdout);
            for (uint32_t j = 0; j < consensus.size(); ++j) {
                fprintf(stdout, ",%u", coverage[i * consensus.size() + j]);
            }
            fputc('\n', stdout);
        }
        // coverage for deletion
        fputc('-', stdout);
        for (uint32_t j = 0; j < consensus.size(); ++j) {
            fprintf(stdout, ",%u", coverage[graph.num_codes() * consensus.size() + j]);
        }
        fputc('\n', stdout);
    }
    else {
        std::string consensus = graph.GenerateConsensus();
        fprintf(stdout, "%s\n%s\n", name.c_str(), consensus.c_str());
    }

    // generate the multiple sequence alignemnt if desired
    if (opt->msa) {
        std::vector<std::string> msa = graph.GenerateMultipleSequenceAlignment(true);
        // remove the extra consensus if added
        if (opt->pairwise_msa) msa.erase(msa.begin());
        // left align if necessary
        if (opt->left_align) left_align_msa(msa);
        for (const auto& it: msa) {
            std::string sequence = it;
            fprintf(stdout, "%s\n", sequence.c_str());
        }
    }
}

void version()
{
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
}

void help(consensus_opt_t *opt)
{
    fprintf(stderr, "Usage: callerpp [options] <sequences>\n\n");
    fprintf(stderr, "A tool for consensus calling and multiple sequence alignment.\n\n");
    fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
    fprintf(stderr, "Input:\n");
    fprintf(stderr, "       The input will be read from standard input. Each batch of\n");
    fprintf(stderr, "       input sequences should start with a FASTA header line ('>').\n");
    fprintf(stderr, "       Each sequence in the batch should be on its own line thereafter.\n\n");
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "       The output will be written to standard output. The input\n");
    fprintf(stderr, "       header line will be written first, followed by the consensus\n");
    fprintf(stderr, "       call, and the multiple sequence alignment if specified.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "       -i, --input     FILE  Read from this input file, stdin otherwise [%s]\n", opt->input.empty() ? "stdin" : opt->input.c_str());
    fprintf(stderr, "       -A, --match     INT   The score for a sequence match [%d]\n", opt->match);
    fprintf(stderr, "       -B, --mismatch  INT   The penalty for a sequence mismatch [%d]\n", opt->mismatch);
    fprintf(stderr, "       -O, --gap       INT   The penalty for a gap [%d]\n", opt->gap);
    fprintf(stderr, "       -a, --algorithm INT   The type of alignment to perform [%d]\n", opt->algorithm);
    fprintf(stderr, "                             0 - local (Smith-Waterman\n");
    fprintf(stderr, "                             1 - global (Needleman-Wunsch)\n");
    fprintf(stderr, "                             2 - semi-global (glocal)\n");
    fprintf(stderr, "       -r, --resort INT      Resort the input sequences prior to POA/MSA [%d]\n", opt->resort);
    fprintf(stderr, "                             0 - do not sort\n");
    fprintf(stderr, "                             1 - by length sequence (shortest first)\n");
    fprintf(stderr, "                             2 - by length sequence (longest first)\n");
    fprintf(stderr, "       -c, --coverage        Output the per-base coverage for the consensus [%s]\n", opt->coverage ? "true" : "false");
    fprintf(stderr, "       -m, --msa             Output multiple sequence alignment [%s]\n", opt->msa ? "true" : "false");
    fprintf(stderr, "       -l, --left-align      Left align the sequences in the multiple sequence alignment [%s]\n", opt->left_align ? "true" : "false");
    fprintf(stderr, "       -p, --pairwise-msa    Re-compute the MSA by adding in the consensus first [%s]\n", opt->pairwise_msa ? "true" : "false");
    fprintf(stderr, "       -h, --help            Prints out the help\n");
    fprintf(stderr, "       -v, --version         Prints out the version\n");
}

int main(int argc, char** argv) {
    int c;
    consensus_opt_t *opt = consensus_opt_init();
    std::vector<std::string> sequences = {};
    std::string name;
    std::unique_ptr<spoa::AlignmentEngine> alignment_engine;
    std::ifstream in;
    std::istream *stream = &std::cin;

    while ((c = getopt_long(argc, argv, "i:A:B:O:a:r:cmlphv", options, nullptr)) != -1) {
        if ('i' == c)      opt->input        = optarg;
        else if ('A' == c) opt->match        = atoi(optarg);
        else if ('B' == c) opt->mismatch     = atoi(optarg);
        else if ('O' == c) opt->gap          = atoi(optarg);
        else if ('a' == c) opt->algorithm    = atoi(optarg);
        else if ('r' == c) opt->resort       = atoi(optarg);
        else if ('c' == c) opt->coverage     = true;
        else if ('m' == c) opt->msa          = true;
        else if ('l' == c) opt->left_align   = true;
        else if ('p' == c) opt->pairwise_msa = true;
        else if ('v' == c) {
            version();
            return 0;
        }
        else {
            help(opt);
            return c == 'h' ? 0 : -1;
        }
    }
    if (optind != argc) {
        help(opt);
        fprintf(stderr, "Error: extra arguments found:");
        while (optind < argc) {
            fprintf(stderr, " %s", argv[optind]);
            optind++;
        }
        fputc('\n', stderr);
        return -1;
    }

    if (opt->input.empty()) {
        fprintf(stderr, "Reading from standard input...\n");
    }
    else {
        in.open(opt->input.c_str(), std::ifstream::in);
        stream = &in;
    }

    alignment_engine = spoa::AlignmentEngine::Create(static_cast<spoa::AlignmentType>(opt->algorithm), opt->match, opt->mismatch, opt->gap);
    std::getline(*stream, name);
    if (name.compare(0, 1, ">") != 0) {
        fprintf(stderr, "Expecting a header line, found: %s\n", name.c_str());
        return 1;
    }
    for (std::string line; std::getline(*stream, line);) {
        // An empty line can trigger process
        if (line.empty() || line.compare(0, 1, ">") == 0) {
            if (sequences.empty()) {
                if (!name.empty()) {
                    fprintf(stderr, "No sequences specified for '%s'\n", name.c_str());
                    return -1;
                }
            }
            else {
                process(alignment_engine, name, sequences, opt);
                if (line.empty()) fflush(stdout);

            }
            name = line.empty() ? "" : line;
            sequences.clear();
        } else {
            if (name.empty()) {
                fprintf(stderr, "No name found for '%s'\n", line.c_str());
                return -1;
            }
            sequences.push_back(line);
        }
    }
    if (!sequences.empty()) {
        process(alignment_engine, name, sequences, opt);
    }

    if (!opt->input.empty()) {
        in.close();
    }
    free(opt);

    return 0;
}
