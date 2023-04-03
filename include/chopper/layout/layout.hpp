#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <chopper/prefixes.hpp>

namespace chopper::layout
{

// Currently, the layout is structured by user bin.
struct layout
{
    struct max_bin
    {
        std::vector<size_t> previous_TB_indices{}; // identifies the IBF based on upper levels
        size_t id{};                               // the technical bin id that has the maximum kmer content

        friend constexpr auto operator<=>(max_bin const &, max_bin const &) = default;

        // needs a template (instead of using std::ostream directly) to be able to only include <iosfwd>
        template <typename stream_type>
            requires std::derived_from<stream_type, std::ostream>
        friend stream_type & operator<<(stream_type & stream, max_bin const & object)
        {
            stream << prefix::header << prefix::merged_bin << '_';
            auto it = object.previous_TB_indices.begin();
            auto end = object.previous_TB_indices.end();
            // If not empty, we join with ';'
            if (it != end)
            {
                stream << *it;
                while (++it != end)
                    stream << ';' << *it;
            }
            stream << " max_bin_id:" << object.id;

            return stream;
        }
    };

    struct user_bin
    {
        std::string id{};                          // the id of the user bin. Often its filename.
        std::vector<size_t> previous_TB_indices{}; // previous technical bin indices which refer to merged bin indices.
        size_t number_of_technical_bins{};         // 1 == signle bin, >1 == split_bin
        size_t storage_TB_id{}; // the id of the technical bin that the user bin is actullly stored in

        friend constexpr auto operator<=>(user_bin const &, user_bin const &) = default;

        // needs a template (instead of using std::ostream directly) to be able to only include <iosfwd>
        template <typename stream_type>
            requires std::derived_from<stream_type, std::ostream>
        friend stream_type & operator<<(stream_type & stream, user_bin const & object)
        {
            stream << object.id << '\t';
            for (auto bin : object.previous_TB_indices)
                stream << bin << ';';
            stream << object.storage_TB_id << '\t';
            for ([[maybe_unused]] auto && elem : object.previous_TB_indices) // number of bins per merged level is 1
                stream << "1;";
            stream << object.number_of_technical_bins;

            return stream;
        }
    };

    std::vector<max_bin> max_bins{};
    std::vector<user_bin> user_bins{};
};

} // namespace chopper::layout
