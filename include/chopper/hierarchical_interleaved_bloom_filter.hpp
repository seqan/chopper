#pragma once

#include <seqan3/std/ranges>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

namespace hibf
{

template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;

    // Forward declaration
    template <std::integral value_t>
    class counting_agent_type;

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of a individual Bloom filter.
    using ibf_t = seqan3::interleaved_bloom_filter<data_layout_mode_>;

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
    std::vector<ibf_t> ibf_vector;

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

    /*!\brief Returns a membership_agent to be used for counting.
     */
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent{*this};
    }

    /*!\brief Returns a counting_agent_type to be used for counting.
     */
    template <std::integral value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }

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

    //!\brief Returns the number of managed user bins.
    size_t num_user_bins() const noexcept
    {
        return filenames.size();
    }

    //!\brief Changes the number of managed IBFs.
    void resize_bins(size_t const size)
    {
        bin_to_filename_position.resize(size);
    }

    //!\brief Changes the number of managed user bins.
    void resize_filename(size_t const size)
    {
        filenames.resize(size);
    }

    //!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
    std::vector<int64_t> & bin_at(size_t const idx)
    {
        return bin_to_filename_position[idx];
    }

    //!\brief Returns the filename of the `idx`th user bin.
    std::string & filename_at(size_t const idx)
    {
        return filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return filenames[bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`the IBF.
              An empty string is returned for merged bins.
     */
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

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
    {
        return bin_to_filename_position[ibf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
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

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        archive(filenames);
        archive(bin_to_filename_position);
    }
    //!\endcond
};

template <seqan3::data_layout data_layout_mode> // TODO: value_t as template?
class hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);

        uint16_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() || // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer.emplace_back(current_filename_index);
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default; //!< Defaulted.
    membership_agent(membership_agent const &) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default; //!< Defaulted.
    ~membership_agent() = default; //!< Defaulted.

    /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent(hibf_t const & hibf) :
        hibf_ptr(std::addressof(hibf))
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector<int64_t> result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each user bin for all values in a range.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
    {
        assert(hibf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        bulk_contains_impl(values, 0, threshold);

        std::ranges::sort(result_buffer); // TODO: necessary?

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) && noexcept = delete;
    //!\}
};

template <seqan3::data_layout data_layout_mode>
template <std::integral value_t>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::counting_agent_type
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive bulk counting.
    template <std::ranges::forward_range value_range_t>
    void bulk_count_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<value_t>();
        auto & result = agent.bulk_count(values);

        value_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];
            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_count_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() || // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer[current_filename_index] = sum;
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit counting_agent_type(hibf_t const & hibf) :
        hibf_ptr(std::addressof(hibf)), result_buffer(hibf_ptr->user_bins.num_user_bins())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    seqan3::counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */

    /*!\brief Counts the occurrences in each user bin for all values in a range.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values, size_t const threshold = 1u) & noexcept
    {
        assert(hibf_ptr != nullptr);
        assert(threshold > 0u);
        assert(result_buffer.size() == hibf_ptr->user_bins.num_user_bins());

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, static_cast<value_t>(0u));

        bulk_count_impl(values, 0, threshold);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values, size_t const threshold = 1u) && noexcept = delete;
    //!\}
};

} // namespace hibf

// TODO
// membership agent that returns vector of user bin ids and uses a threshold
// count agent that just counts and returns count for each user bin: !! do not recurse if count is 0 !!
// Maybe an interface with optional threshold for counting?
// It is apparently no efficient to store the counting agents of the IBFs
// Maybe just store the counting agent of the top level ibf?

// NEW
// counting with just checking if 0 is very slow, offer to give a threshold
// membership agent with threshold and without, returns bin indicies
