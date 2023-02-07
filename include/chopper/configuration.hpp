#pragma once

#include <filesystem>

#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <chopper/cereal/path.hpp>

namespace chopper
{

struct configuration
{   /*!\name General Configuration
     * \{
     */
    //!\brief The input file to chopper. Should contain one file path per line.
    std::filesystem::path data_file;

    //!\brief Internal parameter that triggers some verbose debug output.
    bool debug{false};
    //!\}

    /*!\name Configuration of size estimates (chopper::count)
     * \{
     */
    //!\brief The output prefix for `count_filename` and `sketch_directory`.
    std::string output_prefix{"chopper_sketch"};

    //!\brief Set internally by chopper::detail::apply_prefix.
    std::filesystem::path count_filename{};

    //!\brief Set internally by chopper::detail::apply_prefix.
    std::filesystem::path sketch_directory{};

    //!\brief Size estimates used to be able to be clustered by additional info in the file. OUTDATED.
    size_t column_index_to_cluster{1u};

    //!\brief The kmer size to hash the input sequences before computing a HyperLogLog sketch from them.
    uint8_t k{19};

    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};

    //!\brief Do not write the sketches into a dedicated directory.
    bool disable_sketch_output{false};

    //!\brief Whether the input files are precomputed files (.minimizer) instead of sequence files.
    bool precomputed_files{false};
    //!\}

    /*!\name General Configuration
     * \{
     */
    //!\brief Always has to be the same as output prefix.
    std::string input_prefix{};

    //!\brief The name of the layout file to write.
    std::filesystem::path output_filename{"binning.out"};

    //!\brief The maximum number of technical bins on each IBF in the HIBF.
    uint16_t tmax{64};

    //!\brief Analog to `column_index_to_cluster`. OUTDATED
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

    //!\brief Whether the program should determine the best number of IBF bins by doing multiple binning runs.
    bool determine_best_tmax{false};

    //!\brief Whether the programm should compute all binnings up to the given t_max.
    bool force_all_binnings{false};

    //!\brief Whether to print verbose output when computing the statistics when computing the layout.
    bool output_verbose_statistics{false};

    //!\brief The percentage of empty bins sampled during layout computation.
    double update_ubs{0};

    //!\brief Whether update operations of type sequence insertions are expected.
    bool update_seqs{false};
    //!\}

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t const version{2};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(data_file));
        archive(CEREAL_NVP(debug));
        archive(CEREAL_NVP(count_filename));
        archive(CEREAL_NVP(sketch_directory));
        // archive(CEREAL_NVP(column_index_to_cluster));
        archive(CEREAL_NVP(k));
        archive(CEREAL_NVP(sketch_bits));
        archive(CEREAL_NVP(disable_sketch_output));
        archive(CEREAL_NVP(precomputed_files));

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

} // namespace chopper
