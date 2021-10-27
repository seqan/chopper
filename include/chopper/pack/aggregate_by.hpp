#pragma once

#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/utility/views/to.hpp>

#include <chopper/pack/data_store.hpp>

namespace chopper::pack
{

inline void sort_by(data_store & data, int8_t column_to_sort_by)
{
    // Note: We may only sort by extra information

    // generate permutation of indices sorted in descinding order by the sequence lengths
    auto permutation = std::views::iota(0u, data.filenames.size()) | seqan3::views::to<std::vector>;
    assert(permutation.size() == data.filenames.size());

    auto const & info = data.extra_information;
    auto compare = [&info,column_to_sort_by] (auto const l, auto const r)
        { return info[l][column_to_sort_by] < info[r][column_to_sort_by]; };
    std::sort(permutation.begin(), permutation.end(), compare);

    // apply permutation
    for (size_t i = 0; i < permutation.size(); i++)
    {
        auto current = i;
        while (i != permutation[current])
        {
            auto next = permutation[current];
            std::swap(data.filenames[current], data.filenames[next]);
            std::swap(data.kmer_counts[current], data.kmer_counts[next]);
            std::swap(data.extra_information[current], data.extra_information[next]);
            permutation[current] = current;
            current = next;
        }
        permutation[current] = current;
    }
}

/*!\brief Aggregates data.filenames, dara.kmer_counts and data.extra_information by the specified column index.
 * \param data The data to aggregate.
 * \param column_to_aggregate_by column index referring to the index in data.extra_information (!).
 *
 * \attention The column index to aggregate by should be given relative to the position in data.extra_information.
 *            This differs from the column index the user specifies via the command line because the user input
 *            is a file containing filenames and kmer counts at the first and second column.
 */
inline void aggregate_by(data_store & data, int8_t column_to_aggregate_by)
{
    if (data.filenames.empty())
        return;

    // Note: We may only aggregate by extra information
    assert(data.filenames.size() == data.kmer_counts.size());
    assert(data.filenames.size() == data.extra_information.size());
    assert(column_to_aggregate_by >= 0);
    assert(column_to_aggregate_by < static_cast<int8_t>(data.extra_information[0].size()));

    sort_by(data, column_to_aggregate_by);

    assert(!data.extra_information.empty());
    assert(!data.extra_information[column_to_aggregate_by].empty());

    size_t index{};
    std::string current_info{data.extra_information[index][column_to_aggregate_by]};
    ++index;
    while (index < data.extra_information.size())
    {
        if (current_info == data.extra_information[index][column_to_aggregate_by])
        {
            data.filenames[index - 1] = data.filenames[index - 1] + ";" + data.filenames[index];
            data.kmer_counts[index - 1] = data.kmer_counts[index - 1] + data.kmer_counts[index]; // todo sum is bad
            data.filenames.erase(data.filenames.begin() + index);
            data.kmer_counts.erase(data.kmer_counts.begin() + index);
            data.extra_information.erase(data.extra_information.begin() + index);
        }
        else
        {
            current_info = data.extra_information[index][column_to_aggregate_by];
            ++index;
        }
    }
}

} // namespace chopper::pack
