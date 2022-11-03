#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/count/configuration.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/print_peak_memory_usage.hpp>

namespace chopper::count
{

void initialize_argument_parser(seqan3::argument_parser & parser, chopper::count::configuration & config)
{ //myrthe test 3
    parser.info.author = "Avenja";
    parser.info.short_description = "Count all kmers of each file in a directory.";
    parser.info.version = "1.0.0";

    parser.add_option(config.data_file, 'f', "data_file", "Give me a filename to a seqinfo file.", seqan3::option_spec::required);
    parser.add_option(config.output_filename, 'o', "outfile", "Filename where the counts should be written to.",
                      seqan3::option_spec::required);
    parser.add_option(config.hll_dir, 'd', "hll-dir",
                      "A directory path where the HyperLogLog sketches will be stored in files. "
                      "When this is not given, the sketches will not be stored in files.",
                      seqan3::option_spec::standard,
                      seqan3::output_directory_validator{});
    parser.add_option(config.column_index_to_cluster, 'c', "column-index", "The column index by which to cluster.");
    parser.add_option(config.num_threads, 't', "threads", "Number of threads.");
    parser.add_option(config.k, 'k', "kmer-size", "The kmer size to count minimisers.");
    parser.add_option(config.w, 'w', "window-size", "The window size for minimisers.");
    parser.add_option(config.sketch_bits, 's', "sketch-bits",
                      "The number of bits the HyperLogLog sketch should use to distribute the values into bins.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{5, 32});
    parser.add_flag(config.disable_minimizers, '\0', "disable-minimizers",
                    "Compute pure kmer counts instead of minimizers. Note that selecting -k == -w would not be enough "
                    "because the minimizer hash will still consider the reverse complement and thus differ from a "
                    "pure kmer count.");
    parser.add_flag(config.exclusively_hlls, 'e', "exclusively-hlls",
                    "When this is given, only the Hyperloglog sketches are computed and not exact counts. "
                    "The exact counts are replaced by estimates of the single sketches for every sequence.");
}

int execute(seqan3::argument_parser & parser)
{
    chopper::count::configuration config{};
    initialize_argument_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[CHOPPER COUNT ERROR] " << ext.what() << "\n";
        return -1;
    }

    auto filename_clusters = chopper::count::read_data_file(config);

    chopper::count::count_kmers(filename_clusters, config);

    chopper::print_peak_memory_usage();

    return 0;
}

} // namespace chopper::count
