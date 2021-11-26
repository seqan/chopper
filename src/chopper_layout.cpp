#include <iostream>

#include <seqan3/argument_parser/all.hpp>

#include <chopper/detail_apply_prefix.hpp>
#include <chopper/layout/aggregate_by.hpp>
#include <chopper/layout/configuration.hpp>
#include <chopper/layout/filenames_data_input.hpp>
#include <chopper/layout/hierarchical_binning.hpp>
#include <chopper/layout/ibf_query_cost.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/layout/previous_level.hpp>

namespace chopper::layout
{

void set_up_subparser_layout(seqan3::argument_parser & parser, chopper::layout::configuration & config)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Compute an HIBF layout";

    parser.info.description.emplace_back("Computes an HIBF layout that tries to minimize the disk space consumption of "
                                         "the resulting index. The space is estimated using a k-mer count per user "
                                         "bin which represents the potential denisity in a technical bin in an "
                                         "interleaved Bloom filter.  You can pass the resulting layout to raptor "
                                         "(https://github.com/seqan/raptor) to build the index and "
                                         "conduct queries.");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.input_prefix,
                      '\0', "input-prefix",
                      "Provide the prefix you used for the output prefix in the chopper count --output-prefix option. "
                      "If you have different means of estimating the k-mer counts of your input data, make sure that a "
                      "file [INPUT-PREFIX].count exists. It needs to be tab-separated and consist of two columns: "
                      "\"[filepath] [tab] [weight/count]\".",
                      seqan3::option_spec::required);
    parser.add_list_item("", "Example count file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta     500");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz     600");
    parser.add_list_item("", "```");

    parser.add_option(config.tmax,
                      '\0', "tmax",
                      "Limits the number of technical bins on each level of the HIBF. Choosing a good tmax is not "
                      "trivial. The smaller tmax, the more levels the layout needs to represent the data. This results "
                      "in a higher space consumption of the index. While querying each individual level is cheap, "
                      "querying many levels might also lead to an increased runtime. "
                      "A good tmax is usually the square root of the number of user bins rounded to the next multiple "
                      "of 64. Note that your tmax will be rounded to the next multiple of 64 anyway. "
                      "At the expense of a longer runtime, you can enable the statistic mode that determines the best "
                      "tmax. See the option --determine-best-tmax",
                      seqan3::option_spec::required);

    parser.add_option(config.num_hash_functions,
                      '\0', "num-hash-functions",
                      "The number of hash functions to use when building the HIBF from the resulting layout. "
                      "This parameter is needed to correctly estimate the index size when computing the layout.");

    parser.add_option(config.false_positive_rate,
                      '\0', "false-positive-rate",
                      "The false positive rate you aim for when building the HIBF from the resulting layout. "
                      "This parameter is needed to correctly estimate the index size when computing the layout.");

    parser.add_option(config.output_filename, '\0', "output-filename", "A file name for the resulting layout.");

    using aggregate_by_type = std::remove_cvref_t<decltype(config.aggregate_by_column)>;
    parser.add_option(config.aggregate_by_column,
                      '\0', "aggregate-by-column",
                      "Which column do you want to aggregate your files by? Start counting your columns from 0!",
                      seqan3::option_spec::hidden,
                      seqan3::arithmetic_range_validator{aggregate_by_type{2},
                                                         std::numeric_limits<aggregate_by_type>::max()});

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use. Currently, only merging of sketches is parallelized, so if option "
                      "--rearrange-user-bins is not set, --threads will have no effect.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), std::numeric_limits<size_t>::max()});

    parser.add_subsection("HyperLogLog Sketches:");
    parser.add_line("To improve the layout, you can estimate the sequence similarities using HyperLogLog sketches.");

    parser.add_flag(config.estimate_union,
                    '\0', "estimate-union",
                    "Use sketches to estimate the sequence similarity among a set of user bins. This will improve the "
                    "layout computation as merging user bins that do not increase technical bin sizes will be "
                    "preferred. Attention: Only possible if the directory [INPUT-PREFIX]_sketches is present.");

    parser.add_flag(config.rearrange_user_bins,
                    '\0', "rearrange-user-bins",
                    "As a preprocessing step, rearranging the order of the given user bins based on their sequence "
                    "similarity may lead to favourable small unions and thus a smaller index. "
                    "Attention: Also enables --estimate-union and is only possible if the directory "
                    "[INPUT-PREFIX]_sketches is present.");

    parser.add_subsection("Parameter Tweaking:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.alpha,
                      '\0', "alpha",
                      "The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no "
                      "means of optimizing the runtime for querying such an HIBF. In general, the ratio of merged bins "
                      "and split bins influences the query time because a merged bin always triggers another search on "
                      "a lower level. To influence this ratio, alpha can be used. The higher alpha, the less merged "
                      "bins are chosen in the layout. This improves query times but leads to a bigger index.",
                      seqan3::option_spec::advanced);

    parser.add_option(config.max_rearrangement_ratio,
                      '\0', "max-rearrangement-ratio",
                      "When the option --rearrange-user-bins is set, this option can influence the rearrangement "
                      "algorithm. The algorithm only rearranges the order of user bins in fixed intervals. The higher "
                      "--max-rearrangement-ratio, the larger the intervals. This potentially improves the layout, but "
                      "increases the runtime of the layout algorithm.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0.0, 1.0});

    parser.add_subsection("Special options");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_flag(config.determine_best_tmax,
                    '\0', "determine-best-tmax",
                    "When this flag is set, the program will compute multiple layouts for tmax in "
                    "[64 , 128, 256, ... , tmax] as well as tmax=sqrt(number of user bins). "
                    "The layout algorithm itself only optimizes the space consumption. When determining the best "
                    "layout, we additionally keep track of the average number of queries needed to traverse each "
                    "layout. This query cost is taken into account when determining the best tmax for your data. "
                    "Note that the option --tmax serves as upper bound. Once the layout quality starts dropping, the "
                    "computation is stopped. To run all layout computations, pass the flag --force-all-binnings.");

    parser.add_flag(config.force_all_binnings,
                    '\0', "force-all-binnings",
                    "Forces all layouts up to --tmax to be computed, "
                    "regardless of the layout quality. If the flag --determine-best-tmax is not set, this flag is "
                    "ignored and has no effect.");

    parser.add_flag(config.output_statistics,
                    '\0', "output-statistics",
                    "Enable verbose statistics to be "
                    "printed to std::cout. If the flag --determine-best-tmax is not set, this flag is ignored "
                    "and has no effect.",
                    seqan3::option_spec::advanced);

    parser.add_flag(config.debug,
                    '\0', "debug",
                    "Enables debug output in layout file.",
                    seqan3::option_spec::hidden);
}

void sanity_checks(layout::data_store const & data, chopper::layout::configuration & config)
{
    if (config.rearrange_user_bins)
        config.estimate_union = true;

    if (config.estimate_union &&
        (!std::filesystem::exists(config.sketch_directory) || std::filesystem::is_empty(config.sketch_directory)))
    {
        throw seqan3::argument_parser_error{"The directory " + config.sketch_directory.string() + " must be present "
                                            "and not empty in order to enable --estimate-union or "
                                            "--rearrange-user-bins (created with chopper count)."};
    }

    if (data.filenames.empty())
        throw seqan3::argument_parser_error{seqan3::detail::to_string("The file ", config.count_filename.string(),
                                                                      " appears to be empty.")};

    if (config.aggregate_by_column != -1 && data.extra_information[0].empty())
    {
        throw seqan3::argument_parser_error{"Aggregate Error: You want to aggregate by something but your "
                                            "file does not contain any extra information columns."};
    }

    // note that config.aggregate_by_column cannot be 0 or 1 because the parser check the valid range [2, max]
    assert(config.aggregate_by_column == -1 || config.aggregate_by_column > 1);
    if ((config.aggregate_by_column - 2/*extrainfo starts at 2*/) > static_cast<int>(data.extra_information[0].size()))
    {
        throw seqan3::argument_parser_error{"Aggregate Error: You want to aggregate by a column index that "
                                            "is larger than the number of extra information columns."};
    }
}

size_t determine_best_number_of_technical_bins(chopper::layout::data_store & data, chopper::layout::configuration & config)
{
    std::stringstream * const output_buffer_original = data.output_buffer;
    std::stringstream * const header_buffer_original = data.header_buffer;
    size_t max_hibf_id{};

    // with -determine-best-tmax the algorithm is executed multiple times and result with the minimum
    // expected query costs is written to the output
    std::cout << std::fixed << std::setprecision(2);
    if (!config.output_statistics)
        std::cout << "T_Max\tC_{T_Max}\trelative expected HIBF query cost\n";

    double best_expected_HIBF_query_cost{std::numeric_limits<double>::infinity()};
    size_t best_t_max{};

    size_t const total_kmer_count = std::accumulate(data.kmer_counts.begin(), data.kmer_counts.end(), size_t{});

    std::set<size_t> potential_t_max{};

    for (size_t t_max = 64; t_max <= config.tmax; t_max *= 2)
        potential_t_max.insert(t_max);

    // Additionally, add the t_max that is closest to the sqrt() of the number of
    // user bins, as it is expected to evenly spread bins and may perform well.
    size_t const user_bin_count{std::ranges::size(data.kmer_counts)};
    size_t const sqrt_t_max{next_multiple_of_64(std::ceil(std::sqrt(user_bin_count)))};
    potential_t_max.insert(sqrt_t_max);

    size_t t_max_64_memory{};

    for (size_t const t_max : potential_t_max)
    {
        // reset state for the binning algorithm and save output buffer
        std::stringstream output_buffer_tmp;
        std::stringstream header_buffer_tmp;

        data.output_buffer = &output_buffer_tmp;
        data.header_buffer = &header_buffer_tmp;

        data.previous = chopper::layout::previous_level{};
        config.tmax = t_max;

        chopper::layout::hibf_statistics global_stats{config, data.fp_correction};
        data.stats = &global_stats.top_level_ibf;

        // execute the actual algorithm
        size_t const max_hibf_id_tmp = chopper::layout::hierarchical_binning{data, config}.execute();

        global_stats.finalize();

        double const expected_HIBF_query_cost = global_stats.total_query_cost / total_kmer_count;

        if (!t_max_64_memory)
            t_max_64_memory = global_stats.total_hibf_size_in_byte();

        double const relative_memory_size = global_stats.total_hibf_size_in_byte() /
                                            static_cast<double>(t_max_64_memory);
        double const query_time_memory_usage_prod = expected_HIBF_query_cost * relative_memory_size;

        if (config.output_statistics)
        {
            std::cout << "#T_Max:" << t_max << '\n'
                      << "#C_{T_Max}:" << chopper::layout::ibf_query_cost::interpolated(t_max, config.false_positive_rate) << '\n'
                      << "#relative expected HIBF query time cost (l):" << expected_HIBF_query_cost << '\n' /*relative to a 64 bin IBF*/
                      << "#relative HIBF memory usage (m):" << relative_memory_size << '\n' /*relative to the 64 T_Max HIBF*/
                      << "#l*m:" << query_time_memory_usage_prod << '\n';
        }
        else
        {
            std::cout << t_max << '\t'
                      << chopper::layout::ibf_query_cost::interpolated(t_max, config.false_positive_rate) << '\t'
                      << expected_HIBF_query_cost << '\n';
        }


        if (config.output_statistics)
            global_stats.print_summary();

        // Use result if better than previous one.
        if (expected_HIBF_query_cost < best_expected_HIBF_query_cost)
        {
            *output_buffer_original = std::move(output_buffer_tmp);
            *header_buffer_original = std::move(header_buffer_tmp);
            max_hibf_id = max_hibf_id_tmp;
            best_t_max = t_max;
            best_expected_HIBF_query_cost = expected_HIBF_query_cost;
        }
        else if (!config.force_all_binnings)
        {
            break;
        }
    }

    std::cout << "#Best t_max (regarding expected query runtime):" << best_t_max << '\n';
    return max_hibf_id;
}

int execute(seqan3::argument_parser & parser)
{
    chopper::layout::configuration config;
    chopper::layout::data_store data;

    set_up_subparser_layout(parser, config);

    try
    {
        parser.parse();

        detail::apply_prefix(config.input_prefix, config.count_filename, config.sketch_directory);

        // Read in the data file containing file paths, kmer counts and additional information.
        chopper::layout::read_filename_data_file(data, config);

        sanity_checks(data, config);
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER LAYOUT ERROR] " << ext.what() << '\n';
        return -1;
    }

    if (config.tmax % 64 != 0)
    {
        config.tmax = chopper::next_multiple_of_64(config.tmax);
        std::cerr << "[CHOPPER LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.tmax << ".\n";
    }

    data.compute_fp_correction(config.false_positive_rate, config.num_hash_functions, config.tmax);

    // If requested, aggregate the data before layouting them
    if (config.aggregate_by_column != -1)
        aggregate_by(data, config.aggregate_by_column - 2/*user index includes first two columns (filename, count)*/);

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.false_positive_rate = config.false_positive_rate;

    size_t max_hibf_id;

    // print some important parameters of the user configuration and data input
    if (config.output_statistics)
    {
        std::cout << "#number of user bins:" << data.filenames.size() << '\n'
                  << "#number of hash functions:" << config.num_hash_functions << '\n'
                  << "#false positive rate:" << config.false_positive_rate << "\n\n";
    }

    if (config.determine_best_tmax)
    {
        max_hibf_id = determine_best_number_of_technical_bins(data, config);
    }
    else
    {
        chopper::layout::hibf_statistics global_stats{config, data.fp_correction};
        data.stats = &global_stats.top_level_ibf;

        max_hibf_id = chopper::layout::hierarchical_binning{data, config}.execute(); // just execute once

        if (config.output_statistics)
            global_stats.print_summary();
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    write_layout_header_to(config, max_hibf_id, header_buffer.str(), fout);
    fout << output_buffer.str();

    return 0;
}

} // namespace chopper::layout
