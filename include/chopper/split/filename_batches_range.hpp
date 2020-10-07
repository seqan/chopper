#pragma once

#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

#include <seqan3/range/views/to.hpp>

#include <chopper/detail_starts_with.hpp>
#include <chopper/split/split_config.hpp>

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
        using value_type = split_config;
        //!\brief Pointer type.
        using pointer = void;
        //!\brief Reference type.
        using reference = split_config const &;
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
        iterator(filename_batches_range & host, split_config config) :
            host{&host}, current_config{std::move(config)}
        {
            if (!parse_next_line())
                while (std::getline(this->host->data_file, this->host->current_line) && !parse_next_line());
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
            current_config.seqfiles.clear();
            while (std::getline(host->data_file, host->current_line) && !parse_next_line());
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
            if (host->current_file_type == file_type::seqfiles_given)
                return current_config.seqfiles.empty();
            else
                return host->data_file.eof();
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

        split_config current_config;

        static constexpr std::string_view merged_bin_prefix{"COLORFUL_MERGED_BIN"};

        bool parse_next_line()
        {
            assert(host->current_file_type != file_type::unknown);

            if (host->current_file_type == file_type::seqfiles_given)
                return true;

            char const * buffer = host->current_line.c_str();
            auto field_start = &buffer[0];
            auto field_end = &buffer[0];
            auto const buffer_end = field_start + host->current_line.size();

            while (field_end != buffer_end && *field_end != '\t') ++field_end;

            std::string const bin_name = std::string(field_start, field_end);

            ++field_end; // skip tab
            field_start = field_end;
            while (field_end != buffer_end && *field_end != '\t') ++field_end;

            std::string filenames = std::string(field_start, field_end);

            // read number of technical bins assigned to these files
            ++field_end; // skip tab
            size_t num_technical_bins;
            auto res = std::from_chars(field_end, buffer_end, num_technical_bins);

            if (num_technical_bins == 1)
                return false;

            for (auto && filename : filenames | std::views::split(';'))
                current_config.seqfiles.push_back((filename | seqan3::views::to<std::string>));

            current_config.bins = num_technical_bins;

            std::string out_filename = host->config.out_path.string() + "_" + bin_name + ".out";
            current_config.out_path = std::filesystem::path{out_filename};

            return true;
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

            std::getline(data_file, current_line); // skip header line
            std::getline(data_file, current_line); // read first line
        }

        return identified_file_type;
    }
};

