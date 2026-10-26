// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>
#include <filesystem>
#include <iosfwd>
#include <istream>
#include <ostream>
#include <string>

#include <cereal/cereal.hpp>

#include <hibf/cereal/path.hpp> // IWYU pragma: keep
#include <hibf/config.hpp>
#include <hibf/misc/timer.hpp>

namespace chopper
{

//!\brief The partitioning approaches of the partitioned HIBF. Selected by `configuration::partitioning_approach`.
enum struct partitioning_scheme : uint8_t
{
    blocked,       // 0
    sorted,        // 1
    folded,        // 2
    weighted_fold, // 3
    similarity,    // 4
    lsh,           // 5
    lsh_sim,       // 6
    unknown        // for error parsing
};

struct configuration
{
    //!\brief Whether to use the fast layout algorithm instead of the default one.
    bool fast_layout{false};

    /*!\name General Configuration
     * \{
     */
    //!\brief The input file to chopper. Should contain one file path per line.
    std::filesystem::path data_file;

    //!\brief Internal parameter that triggers some verbose debug output.
    bool debug{false};

    //!\brief The name of the layout file to write.
    std::filesystem::path output_filename{"layout.txt"};

    //!\brief If specified, layout timings are written to the specified file.
    std::filesystem::path output_timings{};

    //!\brief The kmer size to hash the input sequences before computing a HyperLogLog sketch from them.
    uint8_t k{19};

    //!\brief The window size to compute minimizers before computing a HyperLogLog sketch from them.
    uint8_t window_size{k};

    //!\brief Whether the input files are precomputed files (.minimiser) instead of sequence files.
    bool precomputed_files{false};
    //!\}

    /*!\name Partitioned HIBF configuration
     * \{
     */
    //!\brief The number of partitions for the HIBF index. 1 computes a single HIBF layout.
    size_t number_of_partitions{1};

    //!\brief The partitioning approach. See chopper::partitioning_scheme.
    partitioning_scheme partitioning_approach{};

    /*!\brief Whether `hibf_config.tmax` was given by the user. Not serialised.
     *
     * If set, every partition of a partitioned HIBF is laid out with `hibf_config.tmax`. Otherwise, the tmax of each
     * partition is chosen based on its number of user bins.
     */
    bool tmax_is_set{false};
    //!\}

    /*!\name Configuration of size estimates
     * \{
     */
    //!\brief The name for the output directory when writing sketches to disk.
    std::filesystem::path sketch_directory{};

    //!\brief Do not write the sketches into a dedicated directory.
    bool disable_sketch_output{false};
    //!\}

    /*!\name Statistics configuration
     * \{
     */
    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs.
    bool determine_best_tmax{false};

    //!\brief Whether the programm should compute all binnings up to the given t_max.
    bool force_all_binnings{false};

    //!\brief Whether to print verbose output when computing the statistics when computing the layout.
    bool output_verbose_statistics{false};
    //!\}

    //!\brief The HIBF config which will be used to compute the layout within the HIBF lib.
    seqan::hibf::config hibf_config;

    mutable seqan::hibf::concurrent_timer compute_sketches_timer{};
    mutable seqan::hibf::concurrent_timer union_estimation_timer{};
    mutable seqan::hibf::concurrent_timer rearrangement_timer{};
    mutable seqan::hibf::concurrent_timer dp_algorithm_timer{};
    /*!\brief Fast layout: time spent in LSH clustering (`lsh_in_seconds` in the timing output).
     *
     * Summed over all calls of lsh_distributed_ibf_layout, including concurrent ones, so it can exceed the wall-clock
     * time.
     */
    mutable seqan::hibf::concurrent_timer lsh_algorithm_timer{};
    /*!\brief Fast layout: time spent assigning clusters to technical bins by similarity (`search_best_p_in_seconds` in
     *        the timing output).
     *
     * Summed over all calls of lsh_distributed_ibf_layout, including concurrent ones, so it can exceed the wall-clock
     * time.
     */
    mutable seqan::hibf::concurrent_timer find_best_technical_bin_algorithm_timer{};
    /*!\brief Fast layout: time for distributing the user bins onto the technical bins of the top level
     *        (`top_level_lsh_distribution_timer_in_seconds` in the timing output).
     */
    mutable seqan::hibf::concurrent_timer top_level_lsh_distribution_timer{};
    //!\brief Fast layout: time for laying out all lower levels (`small_layouts_timer_in_seconds` in the timing output).
    mutable seqan::hibf::concurrent_timer small_layouts_timer{};

    void read_from(std::istream & stream);

    void write_to(std::ostream & stream) const;

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        // Version 3 added fast_layout.
        // Version 4 added the partitioned HIBF configuration.
        uint32_t version{4};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(data_file));
        archive(CEREAL_NVP(debug));
        archive(CEREAL_NVP(sketch_directory));
        archive(CEREAL_NVP(k));
        archive(CEREAL_NVP(window_size));
        archive(CEREAL_NVP(disable_sketch_output));
        archive(CEREAL_NVP(precomputed_files));

        archive(CEREAL_NVP(output_filename));
        archive(CEREAL_NVP(determine_best_tmax));
        archive(CEREAL_NVP(force_all_binnings));

        if (version >= 3)
            archive(CEREAL_NVP(fast_layout));

        if (version >= 4)
        {
            archive(CEREAL_NVP(number_of_partitions));
            archive(CEREAL_NVP(partitioning_approach));
        }
    }
};

} // namespace chopper

namespace chopper
{

// Make chopper::partitioning_scheme satisfy sharg::ostreamable. Found by argument-dependent lookup.
inline std::ostream & operator<<(std::ostream & output, chopper::partitioning_scheme const & scheme)
{
    if (scheme == chopper::partitioning_scheme::blocked)
        output << "blocked";
    else if (scheme == chopper::partitioning_scheme::sorted)
        output << "sorted";
    else if (scheme == chopper::partitioning_scheme::folded)
        output << "folded";
    else if (scheme == chopper::partitioning_scheme::weighted_fold)
        output << "weighted_fold";
    else if (scheme == chopper::partitioning_scheme::similarity)
        output << "similarity";
    else if (scheme == chopper::partitioning_scheme::lsh)
        output << "lsh";
    else if (scheme == chopper::partitioning_scheme::lsh_sim)
        output << "lsh_sim";
    else
        output << "unknown";

    return output;
}

// Make chopper::partitioning_scheme satisfy sharg::istreamable. Found by argument-dependent lookup.
inline std::istream & operator>>(std::istream & input, chopper::partitioning_scheme & scheme)
{
    std::string text;
    input >> text;

    if (text == "blocked")
        scheme = chopper::partitioning_scheme::blocked;
    else if (text == "sorted")
        scheme = chopper::partitioning_scheme::sorted;
    else if (text == "folded")
        scheme = chopper::partitioning_scheme::folded;
    else if (text == "weighted_fold")
        scheme = chopper::partitioning_scheme::weighted_fold;
    else if (text == "similarity")
        scheme = chopper::partitioning_scheme::similarity;
    else if (text == "lsh")
        scheme = chopper::partitioning_scheme::lsh;
    else if (text == "lsh_sim")
        scheme = chopper::partitioning_scheme::lsh_sim;
    else
        scheme = chopper::partitioning_scheme::unknown;

    return input;
}

} // namespace chopper
