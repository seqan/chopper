#pragma once

#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_parse_chopper_pack_line.hpp>
#include <chopper/detail_starts_with.hpp>
#include <chopper/split/split_config.hpp>

// implements an input range
class filename_batches_range
{
public:
    enum class file_type
    {
        unknown,
        seqfiles_given,
        data_file
    };

private:
    // declaration
    class iterator
    {

    public:
        /*!\name Associated types
         * \{
         */
        //!\brief Difference type.
        using difference_type = int64_t;
        //!\brief Value type.
        using value_type = batch_config;
        //!\brief Pointer type.
        using pointer = void;
        //!\brief Reference type.
        using reference = batch_config const &;
        //!\brief Iterator category.
        using iterator_category = std::input_iterator_tag;
        //!\}

        /*!\name Construction, destruction and assignment
         * \{
         */
        //!\brief Default construction.
        iterator() = default;
        //!\brief Copy construction.
        iterator(iterator const & rhs) = default;
        //!\brief Move construction.
        iterator(iterator && rhs) = default;
        //!\brief Copy assignment.
        iterator & operator=(iterator const & rhs) = default;
        //!\brief Move assignment.
        iterator & operator=(iterator && rhs) = default;
        //!\brief Destruction.
        ~iterator() = default;

        //!\brief Constructing from the underlying seqan3::single_pass_input_view.
        iterator(filename_batches_range & host, split_config const & config) :
            host{&host}, current_config{config}
        {
            parse_next_line();
        }
        //!\}

        /*!\name Access operations
         * \{
         */
        //!\brief Dereferences the cached iterator.
        reference operator*() const noexcept
        {
            return current_config;
        }
        //!\}

        /*!\name Iterator operations
         * \{
         */
        //!\brief Pre-increment.
        iterator & operator++() noexcept
        {
            at_end = !std::getline(host->data_file, host->current_line) || parse_next_line();
            return *this;
        }

        //!\brief Post-increment.
        auto operator++(int) noexcept
        {
            iterator tmp{*this};
            ++(*this);
            return tmp;
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */
        //!\brief Compares for equality with sentinel.
        bool operator==(std::default_sentinel_t const & s) const noexcept
        {
            return at_end;
        }

        //!\copydoc operator==
        friend bool operator==(std::default_sentinel_t const & s, iterator const & rhs) noexcept
        {
            return rhs == s;
        }

        //!\brief Compares for inequality with sentinel.
        bool operator!=(std::default_sentinel_t const & rhs) const noexcept
        {
            return !(*this == rhs);
        }

        //!\copydoc operator!=
        friend bool operator!=(std::default_sentinel_t const & s, iterator const & rhs) noexcept
        {
            return rhs != s;
        }
        //!\}

    protected:
        filename_batches_range * host{nullptr};

        batch_config current_config;

        bool at_end{false};

        bool parse_next_line()
        {
            assert(host->current_file_type != file_type::unknown);

            current_config = batch_config{host->config}; // reset

            if (host->current_file_type == file_type::seqfiles_given)
            {
                current_config.seqfiles = host->config.seqfiles;
                current_config.bin_indices.push_back(0);
                current_config.bins = host->config.bins;
                return true; // end reached
            }

            assert(!host->current_line.empty());
            auto const pack_record = parse_chopper_pack_line(host->current_line);

            current_config.seqfiles = std::move(pack_record.filenames);
            current_config.bin_indices = pack_record.bin_indices;
            current_config.bins = pack_record.number_of_bins.back();

            if (current_config.bin_indices.size() > 1)
            {
                current_config.merged_bin = true;
            }
            else
            {
                current_config.merged_bin = false;
            }

            return host->data_file.eof(); // end not reached yet
        }
    };

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    filename_batches_range() = default; //!< Defaulted.
    filename_batches_range(filename_batches_range const &) = default; //!< Defaulted.
    filename_batches_range(filename_batches_range &&) = default; //!< Defaulted.
    filename_batches_range & operator=(filename_batches_range const &) = default; //!< Defaulted.
    filename_batches_range & operator=(filename_batches_range &&) = default; //!< Defaulted.
    ~filename_batches_range() = default; //!< Defaulted.

    //!\brief Construct from config file
    filename_batches_range(split_config config_) :
        config{std::move(config_)},
        data_file{config.data_filename},
        current_file_type{init()}
    {}
    //!\}

    iterator begin()
    {
        return iterator{*this, config};
    }

    std::default_sentinel_t end()
    {
        return std::default_sentinel_t{};
    }

private:
    split_config config;

    std::ifstream data_file;

    std::string current_line{""};

public:
    file_type const current_file_type{file_type::unknown};

private:

    file_type init()
    {
        file_type identified_file_type{file_type::unknown};

        if (!config.seqfiles.empty())
        {
            identified_file_type = file_type::seqfiles_given;
        }
        else // empty
        {
            identified_file_type = file_type::data_file;

            if (!data_file.good() || !data_file.is_open())
                throw std::logic_error{"Could not open file for reading"};

            while (std::getline(data_file, current_line) && current_line[0] == '#');
        }

        return identified_file_type;
    }
};

