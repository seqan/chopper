#include <algorithm>
#include <limits>

#include <seqan3/argument_parser/all.hpp>

#include <chopper/pack/aggregate_by.hpp>
#include <chopper/pack/filenames_data_input.hpp>
#include <chopper/pack/hierarchical_binning.hpp>
#include <chopper/pack/pack_config.hpp>
#include <chopper/pack/previous_level.hpp>
#include <chopper/union/user_bin_sequence.hpp>
#include <chopper/print_peak_memory_usage.hpp>

#include <robin_hood.h>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, pack_config & config)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    parser.info.description.emplace_back("The `pack` submodule will create a hierarchical binning that minimizes the "
                                         "space consumption of the resulting Interleaved Bloom Filter that you may "
                                         "build with the `build` submodule using the results.");

    parser.add_option(config.data_file, 'f', "filenames",
                      "A tab separated file that contains the filepaths of sequence data you want to analyse.\n"
                      "The first column must contain the paths to sequence files separated by ';'.\n"
                      "The second column must contain the (kmer) count you want you data to be packed into bins. "
                      " See the submodule count for more details on how to add kmer counts to your sequences\n."
                      "All other columns are optional and can be used to aggregate your data (e.g. taxonmic ids).",
                      seqan3::option_spec::required);

    parser.add_option(config.t_max, 'b', "technical-bins",
                      "Into how many technical bins do you want your sequence data to be packed? "
                      "Will be ceiled to the next multiple of 64. Will be an upper bound for the number of bins "
                      "when -determine-num-bins is given.");

    parser.add_option(config.num_hash_functions, 's', "num-hash-functions",
                      "The number of hash functions for the IBFs.");

    parser.add_option(config.fp_rate, 'p', "false-positive-rate",
                      "The desired false positive rate of the IBFs.");

    parser.add_option(config.alpha, 'a', "alpha",
                      "The scaling factor to influence the number of merged bins.");

    parser.add_option(config.max_ratio, 'm', "max-ratio",
                      "[HLL] The maximal cardinality ratio in the clustering intervals.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0.0, 1.0});

    parser.add_option(config.num_threads, 't', "num-threads",
                      "[HLL] The number of threads to use to compute merged HLL sketches.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), std::numeric_limits<size_t>::max()});

    parser.add_option(config.output_filename, 'o', "outfile",
                      "An output file name for the binning results.");

    parser.add_option(config.aggregate_by_column, 'y', "aggregate-by",
                      "Which column do you want to aggregate your files by? Start counting your columns from 1!",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{3, std::numeric_limits<int>::max()});

    parser.add_flag(config.estimate_union, 'u', "estimate-union",
                    "[HLL] Estimate the union of kmer sets to possibly improve the binning.");

    parser.add_option(config.hll_dir, 'd', "hll-dir",
                      "[HLL] If given, the hll sketches are restored from this directory. Required for -u.",
                      seqan3::option_spec::standard,
                      seqan3::input_directory_validator{});

    parser.add_flag(config.rearrange_bins, 'r', "rearrange-bins",
                    "[HLL] Do a rearrangement of the bins which takes into account similarity. Enabling this option "
                    "also enables -u.");

    parser.add_flag(config.determine_num_bins, '\0', "determine-num-bins",
                    "If given, the programm will determine the best number of technical bins by "
                    "doing multiple binning runs. The -b option will then be an upper bound.");

    parser.add_flag(config.debug, '\0', "debug",
                    "Enables debug output in packing file.",
                    seqan3::option_spec::advanced);

    try
    {
        parser.parse();
        if (config.rearrange_bins)
            config.estimate_union = true;
        if (config.estimate_union && !parser.is_option_set("hll-dir") && !parser.is_option_set('d'))
            throw seqan3::argument_parser_error{"An hll dir needs to be provided when enabling -u or -r."};
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER PACK ERROR] " << ext.what() << '\n'; // customize your error message
        return -2;
    }

    return 0;
}

int chopper_pack(seqan3::argument_parser & parser)
{
    pack_config config;

    set_up_and_parse_subparser_split(parser, config);

    // Read in the data file containing file paths, kmer counts and additional information.
    pack_data data;
    read_filename_data_file(data, config);

    config.t_max = next_multiple_of_64(config.t_max);

    // Some sanity checks on the input file and user options.
    if (data.filenames.empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] File Error: you passed an empty file."};

    if (config.aggregate_by_column != -1 && data.extra_information[0].empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] Aggregate Error: You want to aggregate by something but your "
                                 "file does not contain any extra information columns."};
    if (config.aggregate_by_column > static_cast<int>(data.extra_information[0].size()))
        throw std::runtime_error{"[CHOPPER PACK ERROR] Aggregate Error: You want to aggregate by a column index that is "
                                 "larger than the number of extra information columns."};

    // If requested, aggregate the data before packing them
    if (config.aggregate_by_column != -1)
        aggregate_by(data, config.aggregate_by_column);

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    size_t max_hibf_id;

    if (config.determine_num_bins)
    {
        // with -determine-num-bins the algorithm is executed multiple times and result with the minimum
        // expected query costs is written to the output

        // (19,19) mean c_{T_max} (see paper)
        const robin_hood::unordered_map<uint32_t, double> IBF_query_costs = 
        {
            {64, 1.0}, {128, 1.1}, {256, 1.32}, {512, 1.61},
            {1024, 2.69}, {2048, 4.48}, {4096, 7.53}, {8192, 13.65},
            {16384, 23.86}, {32768, 45.66}, {65536, 90.42}
        };

        size_t const total_t_max = config.t_max;
        double best_expected_HIBF_query_cost = std::numeric_limits<double>::max();

        for (size_t t_max = 64; t_max <= total_t_max; t_max *= 2) 
        {
            // reset state for the binning algorithm and save output buffer
            std::stringstream output_buffer_tmp;
            std::stringstream header_buffer_tmp;

            data.output_buffer = &output_buffer_tmp;
            data.header_buffer = &header_buffer_tmp;

            data.previous = previous_level{};
            config.t_max = t_max;

            // A map to keep track of how many k-mers are in split bins per level
            robin_hood::unordered_map<size_t, size_t> kmers_in_split_bins;
            // execute the actual algorithm
            size_t const max_hibf_id_tmp = hierarchical_binning{data, config, kmers_in_split_bins}.execute();

            // calculate the expected number of queries
            size_t weighted_query_nums{};
            size_t total_sum{};
            for (auto && item : kmers_in_split_bins)
            {
                total_sum += item.second;
                weighted_query_nums += (item.first + 1) * item.second;
            }
            
            double const expected_num_queries = static_cast<double>(weighted_query_nums) / total_sum;
            double const expected_HIBF_query_cost = expected_num_queries * IBF_query_costs.at(t_max);

            std::cout << "T_max: " << t_max << " l_{T_max} " << expected_num_queries 
                      << " total " << expected_HIBF_query_cost << '\n';
            
            // check if this is the current best t_max
            if (expected_HIBF_query_cost < best_expected_HIBF_query_cost)
            {
                output_buffer = std::move(output_buffer_tmp);
                header_buffer = std::move(header_buffer_tmp);
                max_hibf_id = max_hibf_id_tmp;
            }
        }
    }
    else 
    {
        // without -determine-num-bins, the binning algorithm is executed just once
        data.output_buffer = &output_buffer;
        data.header_buffer = &header_buffer;

        robin_hood::unordered_map<size_t, size_t> kmers_in_split_bins;
        max_hibf_id = hierarchical_binning{data, config, kmers_in_split_bins}.execute();
    }

    // brief Write the output to the result file.
    std::ofstream fout{config.output_filename};
    fout << "#" << hibf_prefix << " max_bin_id:" << max_hibf_id << '\n';
    fout << header_buffer.str();
    fout << output_buffer.str();

    print_peak_memory_usage();

    return 0;
}
