#pragma once

#include <seqan3/std/ranges>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    class user_bins;

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

    //!\brief The individual interleaved Bloom filters.
    std::vector<seqan3::interleaved_bloom_filter<data_layout_mode_>> ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the id of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

    //!\brief Stores the user bins.
    user_bins user_bins;

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
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond
};


template <seqan3::data_layout data_layout_mode>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::user_bins
{
private:
    //!\brief Containes all filenames.
    std::vector<std::string> filenames;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `bin_to_filename_position[i][b]`.
     * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
     * Otherwise, the returned value `j` can be used to access the corresponding filename `filenames[j]`.
     */
    std::vector<std::vector<int64_t>> bin_to_filename_position{};

public:

    void resize_bins(size_t const size)
    {
        bin_to_filename_position.resize(size);
    }

    void resize_filename(size_t const size)
    {
        filenames.resize(size);
    }

    std::vector<int64_t> & bin_at(size_t const idx)
    {
        return bin_to_filename_position[idx];
    }

    std::string & filename_at(size_t const idx)
    {
        return filenames[idx];
    }

    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return filenames[bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    auto operator[](size_t const ibf_idx) const
    {
        return bin_to_filename_position[ibf_idx]
               | std::views::transform([this] (int64_t i)
                 {
                    if (i == -1)
                        return std::string{};
                    else
                        return filenames[i];
                 });
    }

    int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
    {
        return bin_to_filename_position[ibf_idx][bin_idx];
    }

    template <typename stream_t>
    void write_filenames(stream_t & out_stream) const
    {
        size_t position{};
        std::string line{};
        for (auto const & filename : filenames)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
    }

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        archive(filenames);
        archive(bin_to_filename_position);
    }
};
