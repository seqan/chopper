#include <sharg/parser.hpp>

#include <seqan3/core/debug_stream.hpp>

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>

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

        if (config.update_ubs != 0) // insert empty bins
        {
            std::vector<size_t> insertion_indices = empty_bin_indices(config.update_ubs, data->positions.size());
            empty_bins.resize(sketches_.size());

            insert_empty_bins(insertion_indices, config.estimate_union, config.sketch_bits, filenames, sketches, kmer_counts, empty_bins);
            data.insert_empty_bins(insertion_indices); // filenames and kmer counts are already updated in the line before

            empty_bin_cum_sizes.resize(sketches_.size());

            // create empty_bin_cum_sizes with cumulative sizes
            if (empty_bins[0])
                empty_bin_cum_sizes[0] = kmer_counts.at(0);

            for (size_t j = 1; j < sketches.size(); ++j)
            {
                if (empty_bins[j])
                    empty_bin_cum_sizes[j] = kmer_counts.at(j);
                empty_bin_cum_sizes[j] += empty_bin_cum_sizes[j - 1];
            }
        }

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
