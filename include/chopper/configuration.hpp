#pragma once

#include <filesystem>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <chopper/cereal/path.hpp>

namespace chopper
{

struct configuration
{
    // count configuration
    std::filesystem::path data_file;
    std::string output_prefix{"chopper_sketch"}; // given by the user
    std::filesystem::path count_filename{};   // set internally
    std::filesystem::path sketch_directory{}; // set internally
    size_t column_index_to_cluster{1u};
    uint8_t k{19};
    uint8_t sketch_bits{12};
    bool disable_sketch_output{false};
    //!\brief Whether the input files are precomputed files instead of sequence files.
    bool precomputed_files{false};

    // layout configuration
    std::string input_prefix{}; // provided by user
    std::filesystem::path output_filename{"binning.out"};
    uint16_t tmax{64};
    int8_t aggregate_by_column{-1};
    //!\brief The number of hash functions for the IBFs.
    size_t num_hash_functions{2};
    //!\brief The desired false positive rate of the IBFs.
    double false_positive_rate{0.05};
    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double alpha{1.2};
    //!\brief The maximal cardinality ratio in the clustering intervals.
    double max_rearrangement_ratio{0.5};
    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t threads{1u};
    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool estimate_union{false};
    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool rearrange_user_bins{false};
    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs
    bool determine_best_tmax{false};
    //!\brief Whether the programm should compute all binnings up to the given t_max
    bool force_all_binnings{false};
    bool output_verbose_statistics{false};
    bool debug{false};

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t const version{1};

        archive(CEREAL_NVP(version));
        archive(CEREAL_NVP(input_prefix));
        archive(CEREAL_NVP(count_filename));
        archive(CEREAL_NVP(sketch_directory));
        archive(CEREAL_NVP(output_filename));
        archive(CEREAL_NVP(tmax));
        // archive(CEREAL_NVP(aggregate_by_column));
        archive(CEREAL_NVP(num_hash_functions));
        archive(CEREAL_NVP(false_positive_rate));
        archive(CEREAL_NVP(alpha));
        archive(CEREAL_NVP(max_rearrangement_ratio));
        archive(CEREAL_NVP(threads));
        archive(CEREAL_NVP(estimate_union));
        archive(CEREAL_NVP(rearrange_user_bins));
        archive(CEREAL_NVP(determine_best_tmax));
        archive(CEREAL_NVP(force_all_binnings));
    }
};

} // namespace chopper::layout
