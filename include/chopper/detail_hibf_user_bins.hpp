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

#include <chopper/search/sync_out.hpp>

struct hibf_user_bins
{
private:
    std::vector<std::string> filenames;

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

    void write_filenames(sync_out & out_stream) const
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

// CEREAL_CLASS_VERSION(hibf_user_bins, 1);
