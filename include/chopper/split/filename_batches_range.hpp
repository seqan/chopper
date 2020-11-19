#pragma once

#include <chopper/detail_starts_with.hpp>
#include <chopper/detail_parse_binning_line.hpp>
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
        static constexpr size_t merged_bin_prefix_length{merged_bin_prefix.size()};

        bool parse_next_line()
        {
            assert(host->current_file_type != file_type::unknown);

            if (host->current_file_type == file_type::seqfiles_given)
                return true;

            current_config = host->config; // reset

            auto const bin_data = parse_binning_line(host->current_line);

            if (bin_data.bins == 1)
                return false;

            current_config.seqfiles = std::move(bin_data.filenames);
            current_config.bins = bin_data.bins;

            std::string out_filename;

            if (starts_with(bin_data.bin_name, merged_bin_prefix))
            {
                auto const id_end = std::find(bin_data.bin_name.begin() + merged_bin_prefix_length + 1,
                                              bin_data.bin_name.end(),
                                              '_');
                std::string const merged_bin_id{&bin_data.bin_name[0],
                                                static_cast<std::string_view::size_type>(id_end - bin_data.bin_name.begin())};

                auto it = host->colourful_bin_offsets.find(merged_bin_id);

                if (it == host->colourful_bin_offsets.end()) // merged bin has not been seen yet
                    it = host->colourful_bin_offsets.emplace(merged_bin_id, 0u).first;

                current_config.bin_index_offset = it->second;
                it->second += bin_data.bins;

                out_filename = host->config.out_path.string() + std::string{merged_bin_id} + ".out";
            }
            else
            {
                out_filename = host->config.out_path.string() + bin_data.bin_name + ".out";
            }

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

    std::unordered_map<std::string, size_t> colourful_bin_offsets{};

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

