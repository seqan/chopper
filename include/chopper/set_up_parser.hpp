#pragma once

#include <sharg/parser.hpp>

#include <chopper/configuration.hpp>

inline void set_up_parser(sharg::parser & parser, chopper::configuration & config)
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

    parser.add_option(
            config.update_ubs,
            sharg::config{
                    .short_id = '\0',
                    .long_id = "update-UBs",
                    .description =
                    "Provide a percentage of how many extra UBs need to be sampled",
            }); //myrthe

    parser.add_flag(
            config.update_seqs,
            sharg::config{
                    .short_id = '\0',
                    .long_id = "update-sequences",
                    .description =
                    "", //todo
                    }); //myrthe

    parser.add_option(
        config.data_file,
        sharg::config{
            .short_id = '\0',
            .long_id = "input-file",
            .description =
                "The input must be a file containing paths to sequence data you wish to estimate; one filepath "
                "per line. If your file contains auxiliary information (e.g. species IDs), your file must be tab-"
                "separated.",
            .required = true});
    parser.add_list_item("", "Example file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz");
    parser.add_list_item("", "```");

    // clustering (==aggregating is not selectable at the moment)
    parser.add_option(
        config.column_index_to_cluster,
        sharg::config{.short_id = '\0',
                      .long_id = "column-index",
                      .description = "The column index by which to cluster. Clustering is not supported at the moment",
                      .hidden = true});

    parser.add_option(
        config.k,
        sharg::config{
            .short_id = '\0',
            .long_id = "kmer-size",
            .description =
                "The k-mer size influences the size estimates of the input. "
                "Choosing a k-mer size that is too small for "
                "your data will result in files appearing more similar than they really are. Likewise, a large "
                "k-mer size might miss out on certain similarities. For DNA sequences, a k-mer size between "
                "[16,32] has proven to work well."});

    parser.add_option(
        config.tmax,
        sharg::config{
            .short_id = '\0',
            .long_id = "tmax",
            .description =
                "Limits the number of technical bins on each level of the HIBF. Choosing a good tmax is not "
                "trivial. The smaller tmax, the more levels the layout needs to represent the data. This results "
                "in a higher space consumption of the index. While querying each individual level is cheap, "
                "querying many levels might also lead to an increased runtime. "
                "A good tmax is usually the square root of the number of user bins rounded to the next multiple "
                "of 64. Note that your tmax will be rounded to the next multiple of 64 anyway. "
                "At the expense of a longer runtime, you can enable the statistic mode that determines the best "
                "tmax. See the option --determine-best-tmax",
            .required = true});

    parser.add_option(
        config.num_hash_functions,
        sharg::config{.short_id = '\0',
                      .long_id = "num-hash-functions",
                      .description =
                          "The number of hash functions to use when building the HIBF from the resulting layout. "
                          "This parameter is needed to correctly estimate the index size when computing the layout."});

    parser.add_option(
        config.false_positive_rate,
        sharg::config{.short_id = '\0',
                      .long_id = "false-positive-rate",
                      .description =
                          "The false positive rate you aim for when building the HIBF from the resulting layout. "
                          "This parameter is needed to correctly estimate the index size when computing the layout."});

    parser.add_option(config.output_filename,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output-filename",
                                    .description = "A file name for the resulting layout."});

    // using aggregate_by_type = std::remove_cvref_t<decltype(config.aggregate_by_column)>;
    // parser.add_option(config.aggregate_by_column,
    //                   '\0', "aggregate-by-column",
    //                   "Which column do you want to aggregate your files by? Start counting your columns from 0!",
    //                   seqan3::option_spec::hidden,
    //                   seqan3::arithmetic_range_validator{aggregate_by_type{2},
    //                                                      std::numeric_limits<aggregate_by_type>::max()});

    parser.add_option(
        config.threads,
        sharg::config{
            .short_id = '\0',
            .long_id = "threads",
            .description =
                "The number of threads to use. Currently, only merging of sketches is parallelized, so if option "
                "--rearrange-user-bins is not set, --threads will have no effect.",
            .validator =
                sharg::arithmetic_range_validator{static_cast<size_t>(1), std::numeric_limits<size_t>::max()}});

    parser.add_subsection("HyperLogLog Sketches:");
    parser.add_line("To improve the layout, you can estimate the sequence similarities using HyperLogLog sketches.");

    parser.add_flag(
        config.estimate_union,
        sharg::config{
            .short_id = '\0',
            .long_id = "estimate-union",
            .description =
                "Use sketches to estimate the sequence similarity among a set of user bins. This will improve the "
                "layout computation as merging user bins that do not increase technical bin sizes will be "
                "preferred. Attention: Only possible if the directory [INPUT-PREFIX]_sketches is present."});

    parser.add_flag(
        config.rearrange_user_bins,
        sharg::config{
            .short_id = '\0',
            .long_id = "rearrange-user-bins",
            .description =
                "As a preprocessing step, rearranging the order of the given user bins based on their sequence "
                "similarity may lead to favourable small unions and thus a smaller index. "
                "Attention: Also enables --estimate-union and is only possible if the directory "
                "[INPUT-PREFIX]_sketches is present."});

    parser.add_subsection("Parameter Tweaking:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(
        config.alpha,
        sharg::config{
            .short_id = '\0',
            .long_id = "alpha",
            .description =
                "The layout algorithm optimizes the space consumption of the resulting HIBF but currently has no "
                "means of optimizing the runtime for querying such an HIBF. In general, the ratio of merged bins "
                "and split bins influences the query time because a merged bin always triggers another search on "
                "a lower level. To influence this ratio, alpha can be used. The higher alpha, the less merged "
                "bins are chosen in the layout. This improves query times but leads to a bigger index.",
            .advanced = true});

    parser.add_option(
        config.max_rearrangement_ratio,
        sharg::config{
            .short_id = '\0',
            .long_id = "max-rearrangement-ratio",
            .description =
                "When the option --rearrange-user-bins is set, this option can influence the rearrangement "
                "algorithm. The algorithm only rearranges the order of user bins in fixed intervals. The higher "
                "--max-rearrangement-ratio, the larger the intervals. This potentially improves the layout, but "
                "increases the runtime of the layout algorithm.",
            .advanced = true,
            .validator = sharg::arithmetic_range_validator{0.0, 1.0}});

    parser.add_option(
        config.sketch_bits,
        sharg::config{.short_id = '\0',
                      .long_id = "sketch-bits",
                      .description =
                          "The number of bits the HyperLogLog sketch should use to distribute the values into bins.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{5, 32}});

    parser.add_subsection("Special options");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_flag(
        config.determine_best_tmax,
        sharg::config{
            .short_id = '\0',
            .long_id = "determine-best-tmax",
            .description =
                "When this flag is set, the program will compute multiple layouts for tmax in "
                "[64 , 128, 256, ... , tmax] as well as tmax=sqrt(number of user bins). "
                "The layout algorithm itself only optimizes the space consumption. When determining the best "
                "layout, we additionally keep track of the average number of queries needed to traverse each "
                "layout. This query cost is taken into account when determining the best tmax for your data. "
                "Note that the option --tmax serves as upper bound. Once the layout quality starts dropping, the "
                "computation is stopped. To run all layout computations, pass the flag --force-all-binnings."});

    parser.add_flag(
        config.force_all_binnings,
        sharg::config{
            .short_id = '\0',
            .long_id = "force-all-binnings",
            .description =
                "Forces all layouts up to --tmax to be computed, "
                "regardless of the layout quality. If the flag --determine-best-tmax is not set, this flag is "
                "ignored and has no effect."});

    parser.add_flag(
        config.output_verbose_statistics,
        sharg::config{.short_id = '\0',
                      .long_id = "output-verbose-statistics",
                      .description =
                          "Enable verbose statistics to be "
                          "printed to std::cout. If the flag --determine-best-tmax is not set, this flag is ignored "
                          "and has no effect.",
                      .hidden = true});

    parser.add_flag(
        config.disable_sketch_output,
        sharg::config{
            .short_id = '\0',
            .long_id = "disable-sketch-output",
            .description =
                "Although the sketches will improve the layout, you might want to disable "
                "writing the sketch files to disk. Doing so will save disk space. However, you cannot use either "
                "--estimate-unions or --rearrange-user-bins in chopper layout without the sketches. Note that "
                "this option does not decrease run time as sketches have to be computed either way."});

    parser.add_flag(config.debug,
                    sharg::config{.short_id = '\0',
                                  .long_id = "debug",
                                  .description = "Enables debug output in layout file.",
                                  .hidden = true});

    parser.add_section("References");
    parser.add_line("[1] Philippe Flajolet, Éric Fusy, Olivier Gandouet, Frédéric Meunier. HyperLogLog: the analysis "
                    "of a near-optimal cardinality estimation algorithm. AofA: Analysis of Algorithms, Jun 2007, Juan "
                    "les Pins, France. pp.137-156. hal-00406166v2, https://doi.org/10.46298/dmtcs.3545");
}
