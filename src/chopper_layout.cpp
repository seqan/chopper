#include <iostream>

#include <sharg/all.hpp>

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

void set_up_subparser_layout(sharg::parser & parser, chopper::layout::configuration & config)
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
                      sharg::config{.long_id = "input-prefix",
                                    .description = "Provide the prefix you used for the output prefix in the chopper "
                                                    "count --output-prefix option. If you have different means of "
                                                    "estimating the k-mer counts of your input data, make sure that a "
                                                    "file [INPUT-PREFIX].count exists. It needs to be tab-separated "
                                                    "and consist of two columns: \"[filepath] [tab] [weight/count]\".",
                                    .required = true});
    parser.add_list_item("", "Example count file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta     500");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz     600");
    parser.add_list_item("", "```");

    parser.add_option(config.tmax,
                    sharg::config{.long_id = "tmax",
                                  .description = "Limits the number of technical bins on each level of the HIBF. "
                                                 "Choosing a good tmax is not trivial. The smaller tmax, the more "
                                                 "levels the layout needs to represent the data. This results in a "
                                                 "higher space consumption of the index. While querying each "
                                                 "individual level is cheap, querying many levels might also lead to "
                                                 "an increased runtime. A good tmax is usually the square root of the "
                                                 "number of user bins rounded to the next multiple of 64. Note that "
                                                 "your tmax will be rounded to the next multiple of 64 anyway. At the "
                                                 "expense of a longer runtime, you can enable the statistic mode that "
                                                 "determines the best tmax. See the option --determine-best-tmax",
                                  .required = true});

    parser.add_option(config.num_hash_functions,
                      sharg::config{.long_id = "num-hash-functions",
                                    .description = "The number of hash functions to use when building the HIBF from "
                                                   "the resulting layout. This parameter is needed to correctly "
                                                   "estimate the index size when computing the layout."});

    parser.add_option(config.false_positive_rate,
                    sharg::config{.long_id = "false-positive-rate",
                                  .description = "The false positive rate you aim for when building the HIBF from the "
                                                 "resulting layout. This parameter is needed to correctly estimate the "
                                                 "index size when computing the layout."});

    parser.add_option(config.output_filename,
                      sharg::config{.long_id = "output-filename",
                                    .description = "A file name for the resulting layout."});

    using aggregate_by_type = std::remove_cvref_t<decltype(config.aggregate_by_column)>;
    auto aggregate_max = std::numeric_limits<aggregate_by_type>::max();
    parser.add_option(config.aggregate_by_column,
                      sharg::config{.long_id = "aggregate-by-column",
                                    .description = "Which column do you want to aggregate your files by? Start counting "
                                                   "your columns from 0!",
                                    .hidden = true,
                                    .validator = sharg::arithmetic_range_validator{aggregate_by_type{2}, aggregate_max}});

    parser.add_option(config.threads,
                      sharg::config{.long_id = "threads",
                                    .description = "The number of threads to use. Currently, only merging of sketches "
                                                   "is parallelized, so if option --rearrange-user-bins is not set, "
                                                   "--threads will have no effect.",
                                    .validator = sharg::arithmetic_range_validator{static_cast<size_t>(1),
                                                                                   std::numeric_limits<size_t>::max()}});

    parser.add_subsection("HyperLogLog Sketches:");
    parser.add_line("To improve the layout, you can estimate the sequence similarities using HyperLogLog sketches.");

    parser.add_flag(config.estimate_union,
                    sharg::config{.long_id = "estimate-union",
                                  .description = "Use sketches to estimate the sequence similarity among a set of user "
                                                 "bins. This will improve the layout computation as merging user bins "
                                                 "that do not increase technical bin sizes will be preferred. "
                                                 "Attention: Only possible if the directory [INPUT-PREFIX]_sketches "
                                                 "is present."});

    parser.add_flag(config.rearrange_user_bins,
                    sharg::config{.long_id = "rearrange-user-bins",
                                  .description = "As a preprocessing step, rearranging the order of the given user "
                                                 "bins based on their sequence similarity may lead to favourable small "
                                                 "unions and thus a smaller index. Attention: Also enables "
                                                 "--estimate-union and is only possible if the directory "
                                                 "[INPUT-PREFIX]_sketches is present."});

    parser.add_subsection("Parameter Tweaking:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.alpha,
                      sharg::config{.long_id = "alpha",
                                    .description = "The layout algorithm optimizes the space consumption of the "
                                                   "resulting HIBF but currently has no means of optimizing the "
                                                   "runtime for querying such an HIBF. In general, the ratio of merged "
                                                   "bins and split bins influences the query time because a merged bin "
                                                   "always triggers another search on a lower level. To influence this "
                                                   "ratio, alpha can be used. The higher alpha, the less merged bins "
                                                   "are chosen in the layout. This improves query times but leads to "
                                                   "a bigger index.",
                                    .advanced = true});

    parser.add_option(config.max_rearrangement_ratio,
                    sharg::config{.long_id = "max-rearrangement-ratio",
                                  .description = "When the option --rearrange-user-bins is set, this option can "
                                                 "influence the rearrangement algorithm. The algorithm only rearranges "
                                                 "the order of user bins in fixed intervals. The higher "
                                                 "--max-rearrangement-ratio, the larger the intervals. This "
                                                 "potentially improves the layout, but increases the runtime of the "
                                                 "layout algorithm.",
                                  .validator = sharg::arithmetic_range_validator{0.0, 1.0}});

    parser.add_subsection("Special options");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_flag(config.determine_best_tmax,
                    sharg::config{.long_id = "determine-best-tmax",
                                  .description = "When this flag is set, the program will compute multiple layouts for "
                                                 "tmax in [64 , 128, 256, ... , tmax] as well as tmax=sqrt(number of "
                                                 "user bins). The layout algorithm itself only optimizes the space "
                                                 "consumption. When determining the best layout, we additionally keep "
                                                 "track of the average number of queries needed to traverse each "
                                                 "layout. This query cost is taken into account when determining the "
                                                 "best tmax for your data. Note that the option --tmax serves as upper "
                                                 "bound. Once the layout quality starts dropping, the computation is "
                                                 "stopped. To run all layout computations, pass the flag "
                                                 "--force-all-binnings."});

    parser.add_flag(config.force_all_binnings,
                    sharg::config{.long_id = "force-all-binnings",
                                  .description = "Forces all layouts up to --tmax to be computed, regardless of the "
                                                 "layout quality. If the flag --determine-best-tmax is not set, this "
                                                 "flag is ignored and has no effect."});

    parser.add_flag(config.output_verbose_statistics,
                    sharg::config{.long_id = "output-verbose-statistics",
                                  .description = "Enable verbose statistics to be printed to std::cout. If the flag "
                                                 "--determine-best-tmax is not set, this flag is ignored and has no "
                                                 "effect.",
                                  .hidden = true});

    parser.add_flag(config.debug,
                    sharg::config{.long_id = "debug",
                                  .description = "Enables debug output in layout file.",
                                  .hidden = true});
}

void sanity_checks(layout::data_store const & data, chopper::layout::configuration & config)
{
    if (config.rearrange_user_bins)
        config.estimate_union = true;

    if (config.estimate_union &&
        (!std::filesystem::exists(config.sketch_directory) || std::filesystem::is_empty(config.sketch_directory)))
    {
        throw sharg::parser_error{"The directory " + config.sketch_directory.string() + " must be present "
                                  "and not empty in order to enable --estimate-union or "
                                  "--rearrange-user-bins (created with chopper count)."};
    }

    if (data.filenames.empty())
        throw sharg::parser_error{sharg::detail::to_string("The file ", config.count_filename.string(),
                                                           " appears to be empty.")};

    if (config.aggregate_by_column != -1 && data.extra_information[0].empty())
    {
        throw sharg::parser_error{"Aggregate Error: You want to aggregate by something but your "
                                  "file does not contain any extra information columns."};
    }

    // note that config.aggregate_by_column cannot be 0 or 1 because the parser check the valid range [2, max]
    assert(config.aggregate_by_column == -1 || config.aggregate_by_column > 1);
    if ((config.aggregate_by_column - 2/*extrainfo starts at 2*/) > static_cast<int>(data.extra_information[0].size()))
    {
        throw sharg::parser_error{"Aggregate Error: You want to aggregate by a column index that "
                                  "is larger than the number of extra information columns."};
    }
}

size_t determine_best_number_of_technical_bins(chopper::layout::data_store & data, chopper::layout::configuration & config)
{
    std::stringstream * const output_buffer_original = data.output_buffer;
    std::stringstream * const header_buffer_original = data.header_buffer;

    std::set<size_t> potential_t_max = [&] ()
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

    std::cout << "## ### Parameters ###\n"
              << "## number of user bins = " << data.filenames.size() << '\n'
              << "## number of hash functions = " << config.num_hash_functions << '\n'
              << "## false positive rate = " << config.false_positive_rate << '\n';
    hibf_statistics::print_header(config.output_verbose_statistics);

    double best_expected_HIBF_query_cost{std::numeric_limits<double>::infinity()};
    size_t best_t_max{};
    size_t max_hibf_id{};
    size_t t_max_64_memory{};

    for (size_t const t_max : potential_t_max)
    {
        std::stringstream output_buffer_tmp;
        std::stringstream header_buffer_tmp;
        config.tmax = t_max;                               // overwrite tmax
        data.output_buffer = &output_buffer_tmp;           // overwrite buffer
        data.header_buffer = &header_buffer_tmp;           // overwrite buffer
        data.previous = chopper::layout::previous_level{}; // reset previous IBF, s.t. data refers to top level IBF

        chopper::layout::hibf_statistics global_stats{config, data.fp_correction, data.kmer_counts};
        data.stats = &global_stats.top_level_ibf;

        // execute the actual algorithm
        size_t const max_hibf_id_tmp = chopper::layout::hierarchical_binning{data, config}.execute();

        global_stats.finalize();

        global_stats.print_summary(t_max_64_memory, config.output_verbose_statistics);

        // Use result if better than previous one.
        if (global_stats.expected_HIBF_query_cost < best_expected_HIBF_query_cost)
        {
            *output_buffer_original = std::move(output_buffer_tmp);
            *header_buffer_original = std::move(header_buffer_tmp);
            max_hibf_id = max_hibf_id_tmp;
            best_t_max = t_max;
            best_expected_HIBF_query_cost = global_stats.expected_HIBF_query_cost;
        }
        else if (!config.force_all_binnings)
        {
            break;
        }
    }

    std::cout << "# Best t_max (regarding expected query runtime): " << best_t_max << '\n';
    config.tmax = best_t_max;
    return max_hibf_id;
}

int execute(sharg::parser & parser)
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
    catch (sharg::parser_error const & ext) // the user did something wrong
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
            global_stats.print_header();
            global_stats.print_summary(dummy);
        }
    }

    // brief Write the output to the layout file.
    std::ofstream fout{config.output_filename};
    write_layout_header_to(config, max_hibf_id, header_buffer.str(), fout);
    fout << output_buffer.str();

    return 0;
}

} // namespace chopper::layout
