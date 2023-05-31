#include <sharg/parser.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>
#include <chopper/layout/insert_empty_bins.hpp>

int main(int argc, char const * argv[])
{
    sharg::parser parser{"chopper", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    chopper::configuration config;
    set_up_parser(parser, config);

    try
    {
        parser.parse();

        config.disable_sketch_output = !parser.is_option_set("output-sketches-to");
    }
    catch (sharg::parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    int exit_code{};

    chopper::layout::layout hibf_layout{};
    std::vector<std::string> filenames{};
    std::vector<size_t> kmer_counts{};
    std::vector<bool> empty_bins; // A bitvector indicating whether a bin is empty (1) or not (0).
    std::vector<size_t> empty_bin_cum_sizes; // The cumulative k-mer count for the first empty bin until empty bin i
    std::vector<chopper::sketch::hyperloglog> sketches{};

    chopper::sketch::read_data_file(config, filenames);

    try
    {
        exit_code |= chopper::sketch::execute(config, filenames, sketches);
        chopper::sketch::estimate_kmer_counts(sketches, kmer_counts);

        chopper::layout::insert_empty_bins(empty_bins, empty_bin_cum_sizes,
                          kmer_counts, sketches, filenames, config);

        chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                                  .hibf_layout = &hibf_layout,
                                  .kmer_counts = kmer_counts,
                                  .sketches = sketches,
                                  .empty_bins = empty_bins,
                                  .empty_bin_cum_sizes = empty_bin_cum_sizes};

        exit_code |= chopper::layout::execute(config, filenames, store);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
        return -1;
    }
    return exit_code;
}