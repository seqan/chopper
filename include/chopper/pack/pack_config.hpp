#pragma once

#include <seqan3/std/filesystem>

#include <chopper/pack/previous_level.hpp>

struct pack_config
{
    std::filesystem::path data_file;
    std::filesystem::path output_filename{"binning.out"};
    uint16_t t_max{64};
    uint16_t t_min{64};
    int8_t aggregate_by_column{-1};
    //!\brief If given, the hll sketches are dumped to this directory and restored when they already exist.
    std::filesystem::path hll_dir{};
    //!\brief The number of hash functions for the IBFs.
    size_t num_hash_functions{2};
    //!\brief The desired false positive rate of the IBFs.
    double fp_rate{0.05};
    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double alpha{1.2};
    //!\brief The maximal cardinality ratio in the clustering intervals.
    double max_ratio{0.5};
    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t num_threads{1u};
    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool estimate_union{false};
    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool rearrange_bins{false};
    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs
    bool determine_num_bins{false};
    //!\brief Whether the programm should compute all binnings up to the given t_max
    bool force_all_binnings{false};
    bool output_statistics{false};
    bool debug{false};
};
