#include <seqan3/argument_parser/all.hpp>

#include <chopper/pack/aggregate_by.hpp>
#include <chopper/pack/hierarchical_binning.hpp>
#include <chopper/pack/filenames_data_input.hpp>
#include <chopper/pack/pack_config.hpp>
#include <chopper/print_peak_memory_usage.hpp>

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
                      "Will be ceiled to the next multiple of 64.");

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
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;

    // Execute the actual algorithm:
    size_t const max_hibf_id = hierarchical_binning{data, config}.execute();

    // brief Write the output to the result file.
    std::ofstream fout{config.output_filename};
    fout << "#" << hibf_prefix << " max_bin_id:" << max_hibf_id << '\n';
    fout << header_buffer.str();
    fout << output_buffer.str();

    print_peak_memory_usage();

    return 0;
}
