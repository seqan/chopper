#include <iostream>
#include <set>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/aggregate_by.hpp>
#include <chopper/layout/hierarchical_binning.hpp>
#include <chopper/layout/ibf_query_cost.hpp>
#include <chopper/layout/output.hpp>

namespace chopper::layout
{

size_t determine_best_number_of_technical_bins(chopper::data_store & data, chopper::configuration & config)
{
    chopper::layout::layout * original_layout = data.hibf_layout; // cache original layout

    std::set<size_t> potential_t_max = [&]()
    {
        std::set<size_t> result;

        for (size_t t_max = 64; t_max <= config.tmax; t_max *= 2)
            result.insert(t_max);

        // Additionally, add the t_max that is closest to the sqrt() of the number of
        // user bins, as it is expected to evenly spread bins and may perform well.
        size_t const user_bin_count{std::ranges::size(data.kmer_counts)};
        size_t const sqrt_t_max{next_multiple_of_64(std::ceil(std::sqrt(user_bin_count)))};
        result.insert(sqrt_t_max);

        return result;
    }();

    // with -determine-best-tmax the algorithm is executed multiple times and result with the minimum
    // expected query costs are written to the standard output

    std::ofstream file_out{config.output_filename.string() + ".stats"};

    file_out << "## ### Parameters ###\n"
             << "## number of user bins = " << data.filenames.size() << '\n'
             << "## number of hash functions = " << config.num_hash_functions << '\n'
             << "## false positive rate = " << config.false_positive_rate << '\n';
    hibf_statistics::print_header_to(file_out, config.output_verbose_statistics);

    double best_expected_HIBF_query_cost{std::numeric_limits<double>::infinity()};
    size_t best_t_max{};
    size_t max_hibf_id{};
    size_t t_max_64_memory{};

    for (size_t const t_max : potential_t_max)
    {
        chopper::layout::layout tmp_layout{}; // will be rewritten for every tmax
        config.tmax = t_max;                  // overwrite tmax
        data.hibf_layout = &tmp_layout;
        data.previous = chopper::data_store::previous_level{}; // reset previous IBF, s.t. data refers to top level IBF

        chopper::layout::hibf_statistics global_stats{config, data.fp_correction, data.kmer_counts};
        data.stats = &global_stats.top_level_ibf;

        // execute the actual algorithm
        size_t const max_hibf_id_tmp = chopper::layout::hierarchical_binning{data, config}.execute();

        global_stats.finalize();

        global_stats.print_summary_to(t_max_64_memory, file_out, config.output_verbose_statistics);

        // Use result if better than previous one.
        if (global_stats.expected_HIBF_query_cost < best_expected_HIBF_query_cost)
        {
            *original_layout = std::move(tmp_layout);
            max_hibf_id = max_hibf_id_tmp;
            best_t_max = t_max;
            best_expected_HIBF_query_cost = global_stats.expected_HIBF_query_cost;
        }
        else if (!config.force_all_binnings)
        {
            break;
        }
    }

    file_out << "# Best t_max (regarding expected query runtime): " << best_t_max << '\n';
    config.tmax = best_t_max;
    data.hibf_layout = original_layout; // reset original layout
    return max_hibf_id;
}

int execute(chopper::configuration & config, chopper::data_store & data)
{
    if (config.disable_estimate_union)
        config.disable_rearrangement = true;

    if (config.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (size_t number_samples = data.filenames.size(); number_samples >= 1ULL << 32) // sqrt is bigger than uint16_t
            throw std::invalid_argument{"Too many samples. Please set a tmax (see help via `-hh`)."}; // GCOVR_EXCL_LINE
        else
            config.tmax = chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_samples))));
    }
    else if (config.tmax % 64 != 0)
    {
        config.tmax = chopper::next_multiple_of_64(config.tmax);
        std::cerr << "[CHOPPER LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.tmax << ".\n";
    }

    data.compute_fp_correction(config.false_positive_rate, config.num_hash_functions, config.tmax);

    // TODO aggregating is outdated. Has to be reworked when needed
    // If requested, aggregate the data before layouting them
    // if (config.aggregate_by_column != -1)
    //     aggregate_by(data, config.aggregate_by_column - 2/*user index includes first two columns (filename, count)*/);

    size_t max_hibf_id;

    if (config.determine_best_tmax)
    {
        max_hibf_id = determine_best_number_of_technical_bins(data, config);
    }
    else
    {
        chopper::layout::hibf_statistics global_stats{config, data.fp_correction, data.kmer_counts};
        data.stats = &global_stats.top_level_ibf;
        size_t dummy{};

        max_hibf_id = chopper::layout::hierarchical_binning{data, config}.execute(); // just execute once

        if (config.output_verbose_statistics)
        {
            global_stats.print_header_to(std::cout);
            global_stats.print_summary_to(dummy, std::cout);
        }
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    chopper::layout::write_config_to(config, fout);
    chopper::layout::write_layout_header_to(*(data.hibf_layout), max_hibf_id, fout);
    chopper::layout::write_layout_content_to(*(data.hibf_layout), fout);

    return 0;
}

} // namespace chopper::layout
