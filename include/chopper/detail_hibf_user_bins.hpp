#pragma once

#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/to.hpp>

struct hibf_user_bins
{
private:
    std::vector<std::string> filenames;

    std::vector<std::vector<int64_t>> bin_to_filename_position{};

public:

    size_t add_user_bin(std::string const & name)
    {
        filenames.push_back(name);
        return filenames.size() - 1;
    }

    size_t add_user_bin(std::vector<std::string> const & names)
    {
        // concat
        filenames.push_back((names | seqan3::views::join_with(std::string{";"}) | seqan3::views::to<std::string>));
        return filenames.size() - 1;
    }

    void add_user_bin_positions(std::vector<int64_t> const & bin_indices)
    {
        bin_to_filename_position.push_back(bin_indices);
    }

    void prepend_user_bin_positions(std::vector<int64_t> const & bin_indices)
    {
        bin_to_filename_position.insert(bin_to_filename_position.begin(), bin_indices);
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

    int64_t get_position(size_t const ibf_idx, size_t const bin_idx) const
    {
        return bin_to_filename_position[ibf_idx][bin_idx];
    }

    void write_filenames(std::ostream & out_stream) const
    {
        size_t position{};
        for (auto const & filename : filenames)
        {
            out_stream << '#' << position << '\t' << filename << '\n';
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

// CEREAL_CLASS_VERSION(hibf_user_bins, 1);
