#pragma once

#include <string>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/layout.hpp>
#include <chopper/sketch/toolbox.hpp>

namespace chopper
{

struct data_store
{
    /*!\brief Stores information of the previous level of a given IBF.
     *
     * When computing a hierarchical layout, the `data_store` data structure is used with local copies for each IBF
     * within the hierarchical structure of the HIBF. To keep track of the hierarchy, the `previous_level` stores
     * information about the previous level (where the corresponding merged bin is located).
     */
    struct previous_level
    {
        std::vector<size_t> bin_indices{};
        std::string num_of_bins;

        bool empty() const
        {
            assert(bin_indices.empty() == num_of_bins.empty());
            return bin_indices.empty();
        }
    };

    /*!\name References to global instances of the HIBF.
     * \{
     */
    //!\brief The desired maximum false positive rate of the resulting index.
    double const false_positive_rate{};

    //!\brief The layout that is build by layout::hierarchical_binning.
    layout::layout * hibf_layout;

    //!\brief The kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> const & kmer_counts{};

    //!\brief The hyperloglog sketches of all input files to estimate their size and similarities.
    std::vector<sketch::hyperloglog> const & sketches{};

    //!\brief A bitvector indicating whether a bin is empty (1) or not (0).
    std::vector<bool> const & empty_bins;

    //!\brief The cumulative k-mer count for the first empty bin until empty bin i
    std::vector<size_t> const & empty_bin_cum_sizes;
    //!\}

    /*!\name Local Storage one IBF in the HIBF.
     *
     * These member variables change on each IBF of the HIBF s.t. the current IBF can be constructed from
     * the current subset of data. The same data is also used for the top level IBF that holds all the data.
     * \{
     */

    //!\brief The input is sorted and rearranged. To keep track without changing the input we store the positions.
    std::vector<size_t> positions = [this]()
    {
        std::vector<size_t> ps;
        ps.resize(this->kmer_counts.size());
        std::iota(ps.begin(), ps.end(), 0);
        return ps;
    }();

    //!\brief The false positive correction based on fp_rate, num_hash_functions and requested_max_tb.
    std::vector<double> fp_correction{};

    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level previous{};

    //!\brief Matrix of estimates of merged bin cardinalites
    std::vector<uint64_t> union_estimates{};

    bool user_bins_arranged{false};
    //!\}

};

} // namespace chopper
