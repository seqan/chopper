#pragma once

#include <chopper/search/pair_hash.hpp>
#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

inline void compute_kmers(std::vector<size_t> & kmers, seqan3::dna4_vector const & query, search_config const & config)
{
    kmers.clear();
    auto hash_view = query | seqan3::views::kmer_hash(seqan3::ungapped{config.k});
    std::ranges::move(hash_view, std::cpp20::back_inserter(kmers));
}

inline void write_header(search_data const & data, std::ostream & out_stream)
{
    data.user_bins.write_filenames(out_stream);
    out_stream << "#QUERY_NAME\tUSER_BINS\n";
}

inline void write_result(std::unordered_set<std::pair<int32_t, uint32_t>, pair_hash> const & membership_result,
                         std::string const & id,
                         search_data const & data,
                         std::ostream & out_stream)
{
    if (membership_result.empty())
    {
        out_stream << id << '\t' << std::endl;
        return;
    }

    // storing and sorting this is only done for testing purposes.
    // If this turns out to have a significant runtime penalty, it should be removed.
    std::vector<int64_t> result_positions; // TODO allocate this outside of this function
    for (auto const & [ibf_idx, bin_idx] : membership_result)
    {
        assert(data.user_bins.get_position(ibf_idx, bin_idx) > -1);
        result_positions.push_back(data.user_bins.get_position(ibf_idx, bin_idx));
    }
    std::sort(result_positions.begin(), result_positions.end()); // otherwise the result output is not testable

    out_stream << id << '\t';
    for (size_t i = 0; i < result_positions.size() - 1; ++i)
        out_stream << result_positions[i] << ',';
    out_stream << result_positions.back() << std::endl;
}

inline void search(std::unordered_set<std::pair<int32_t, uint32_t>, pair_hash> & membership_result,
                   std::vector<size_t> const & kmers,
                   search_data const & data,
                   search_config const & config,
                   int64_t const ibf_idx)
{
    size_t const kmer_lemma = (kmers.size() > (config.errors + 1) * config.k)
                              ? kmers.size() - (config.errors + 1) * config.k
                              : 0;

    auto counting_agent = data.hibf[ibf_idx].counting_agent<uint16_t>();

    auto const & result = counting_agent.bulk_count(kmers);

    size_t bin = 0;
    size_t sum = 0;

    assert(result.size() > 0);

    while (bin < result.size() - 1)
    {
        sum += result[bin];

        if (data.user_bins.get_position(ibf_idx, bin) < 0 /*merged bin*/ ||
            data.user_bins.get_position(ibf_idx, bin) != data.user_bins.get_position(ibf_idx, bin + 1))
        {
            if (sum >= kmer_lemma)
            {
                int64_t const next_ibf_idx = data.hibf_bin_levels[ibf_idx][bin];
                if (next_ibf_idx != ibf_idx)
                {
                    search(membership_result, kmers, data, config, next_ibf_idx);
                }
                else
                {
                    membership_result.emplace(ibf_idx, bin);
                }

            }
            sum = 0;
        }
        ++bin;
    }

    sum += result[bin];
    if (sum >= kmer_lemma)
    {
        int64_t const next_ibf_idx = data.hibf_bin_levels[ibf_idx][bin];
        if (next_ibf_idx != ibf_idx)
        {
            search(membership_result, kmers, data, config, next_ibf_idx);
        }
        else
        {
            membership_result.emplace(ibf_idx, bin);
        }

    }
}
