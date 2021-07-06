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

    // maps for each ibf in hibf, each bin to the next ibf postition in hibf (if it is a merged bin)
    // or to the same ibf (if it is not a merged bin). You can thereby check if you need to query another
    // lower level IBF.
    std::vector<std::vector<int64_t>> hibf_bin_levels;

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
        archive(hibf_bin_levels);
        archive(user_bins);
    }
    //!\endcond
};
