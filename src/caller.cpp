#include <iostream>
#include <fstream>
#include <getopt.h>
#include "spoa/spoa.hpp"
#include "spoa/alignment_engine.hpp"
using namespace std;

static struct option options[] = {
    {"input", optional_argument, 0, 'i'},
    {"match", required_argument, 0, 'A'},
    {"mismatch", required_argument, 0, 'B'},
    {"gap", required_argument, 0, 'O'},
    {"algorithm", required_argument, 0, 'a'},
    {"msa", no_argument, 0, 'm'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

typedef struct {
	std::string input;
    int8_t match;
    int8_t mismatch;
    int8_t gap;
    int8_t algorithm;
    bool msa;
} consensus_opt_t;

consensus_opt_t *consensus_opt_init()
{
    consensus_opt_t *opt = (consensus_opt_t*)calloc(1, sizeof(consensus_opt_t));
	opt->input     = "";
    opt->match     = 5;
    opt->mismatch  = -4;
    opt->gap       = -8;
    opt->algorithm = 0;
    opt->msa       = false;
    return opt;
}

void process(std::unique_ptr<spoa::AlignmentEngine> &alignment_engine, std::string &name, std::vector<std::string> &sequences, consensus_opt_t *opt) {
    alignment_engine->prealloc(sequences.size(), 4);

    auto graph = spoa::createGraph();

    // add the alignments to the graph
    for (const auto& it: sequences) {
        auto alignment = alignment_engine->align_sequence_with_graph(it, graph);
        graph->add_alignment(alignment, it);
    }

    // call the consensus
    std::string consensus = graph->generate_consensus();
    fprintf(stdout, "%s\n%s\n", name.c_str(), consensus.c_str());

    // generate the multiple sequence alignemnt if desired
    if (opt->msa) {
        std::vector<std::string> msa;
        graph->generate_multiple_sequence_alignment(msa);

        for (const auto& it: msa) {
            fprintf(stdout, "%s\n", it.c_str());
        }
    }
}

void help(consensus_opt_t *opt)
{
    fprintf(stderr, "Usage: consensus [options] <sequences>\n\n");
    fprintf(stderr, "A tool for consensus calling and multiple sequence alignment.\n\n");
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
    fprintf(stderr, "       -m, --msa             Output multiple sequence alignment [%s]\n", opt->msa ? "true" : "false");
    fprintf(stderr, "       -h, --help            Prints out the help\n");
}

int main(int argc, char** argv) {
    char c;
    consensus_opt_t *opt = consensus_opt_init();
    std::vector<std::string> sequences = {};
    std::string name;
    std::unique_ptr<spoa::AlignmentEngine> alignment_engine;
	std::ifstream in;
	std::istream *stream = &std::cin;

    while ((c = getopt_long(argc, argv, "i:A:B:O:a:mh", options, nullptr)) != -1) {
		if ('i' == c) opt->input = optarg;
		else if ('A' == c) opt->match          = atoi(optarg);
        else if ('B' == c) opt->mismatch  = atoi(optarg);
        else if ('O' == c) opt->gap       = atoi(optarg);
        else if ('a' == c) opt->algorithm = atoi(optarg);
        else if ('m' == c) opt->msa       = true;
        else {
            help(opt);
            return -1;
        }
    }

	if (!opt->input.empty()) {
		in.open(opt->input.c_str(), std::ifstream::in);
		stream = &in;
	}

    alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(opt->algorithm), opt->match, opt->mismatch, opt->gap);
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
				if (line.isempty()) fflush(stdout);

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

    free(opt);
	if (!opt->input.empty()) {
		in.close();
	}

    return 0;
}
