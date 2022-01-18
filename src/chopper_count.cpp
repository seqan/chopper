#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <chopper/prefixes.hpp>
#include <chopper/count/configuration.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/print_peak_memory_usage.hpp>

namespace chopper::count
{

void initialize_argument_parser(seqan3::argument_parser & parser, chopper::count::configuration & config)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Estimate the (representative) k-mer count";

    parser.info.description.emplace_back("Estimates the (representative) k-mer count of sequence data by computing "
                                         "HyperLogLog sketches [1] of the input data. "
                                         "Additionally, the HLL sketches are stored in a directory "
                                         "and can be used in computing an HIBF layout with chopper layout.");

    parser.add_subsection("Main Options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.data_file,
                      '\0', "input-file",
                      "The input must be a file containing paths to sequence data you wish to estimate; one filepath "
                      "per line. If your file contains auxiliary information (e.g. species IDs), your file must be tab-"
                      "separated.",
                      seqan3::option_spec::required);
    parser.add_list_item("", "The paths must either lead to sequence files, in which case the sequences are read "
                             "and hash values are computed using the seqan3::views::kmer_hash or have the extension "
                             ".minimizer, in which case it considered a binary encoded file containing hash values.",
                             seqan3::option_spec::advanced);
    parser.add_list_item("", "Example file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz");
    parser.add_list_item("", "```");

    parser.add_option(config.output_prefix,
                      '\0', "output-prefix",
                      "The program creates a `[PREFIX].count` file that contains estimated k-mer counts for each "
                      "sequence file. These counts are used in chopper layout to create an index layout file that can "
                      "be used by the app raptor to build an HIBF. Additionally, a directory [PREFIX]_sketches is "
                      "created that will contain one `.hll` file per sequence file. The sketch files improve the "
                      "layout computation."
                      "Attention: Expects a prefix, thus file extensions are dropped and if given a path (e.g. /tmp/), "
                      "the default prefix " + std::string{prefix::chopper} + " is appended.",
                      seqan3::option_spec::required);

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to be used for parallel processing.");

    parser.add_option(config.column_index_to_cluster,
                      '\0', "column-index",
                      "The column index by which to cluster. Clustering is not supported at the moment",
                      seqan3::option_spec::hidden);

    parser.add_subsection("Parameter tweaking:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.k,
                      '\0', "kmer-size",
                      "The k-mer size influences the estimated counts. Choosing a k-mer size that is too small for "
                      "your data will result in files appearing more similar than they really are. Likewise, a large "
                      "k-mer size might miss out on certain similarities. For DNA sequences, a k-mer size between "
                      "[16,32] has proven to work well.");

    parser.add_option(config.sketch_bits,
                      '\0', "sketch-bits",
                      "The number of bits the HyperLogLog sketch should use to distribute the values into bins.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{5, 32});

    parser.add_flag(config.disable_sketch_output, '\0', "disable-sketch-output",
                    "Although the sketches will improve the layout, you might want to disable "
                    "writing the sketch files to disk. Doing so will save disk space. However, you cannot use either "
                    "--estimate-unions or --rearrange-user-bins in chopper layout without the sketches. Note that "
                    "this option does not decrease run time as sketches have to be computed either way.");

    parser.add_section("References");
    parser.add_line("[1] Philippe Flajolet, Éric Fusy, Olivier Gandouet, Frédéric Meunier. HyperLogLog: the analysis "
                    "of a near-optimal cardinality estimation algorithm. AofA: Analysis of Algorithms, Jun 2007, Juan "
                    "les Pins, France. pp.137-156. hal-00406166v2, https://doi.org/10.46298/dmtcs.3545");
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

    detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);

    chopper::count::count_kmers(filename_clusters, config);

    chopper::print_peak_memory_usage();

    return 0;
}

} // namespace chopper::count
