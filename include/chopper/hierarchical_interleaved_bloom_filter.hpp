#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <chopper/detail_hibf_user_bins.hpp>

template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter() = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter & operator=(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter & operator=(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default; //!< Defaulted.

    //!\}

    std::vector<seqan3::interleaved_bloom_filter<data_layout_mode_>> hibf;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the id of the lower
     * level IBF in hibf.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

    hibf_user_bins user_bins;

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(hibf);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond
};
