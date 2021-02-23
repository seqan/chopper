#pragma once

#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

struct pair_hash
{
    std::size_t operator () (std::pair<int32_t, uint32_t> const & pair) const
    {
        return (static_cast<size_t>(pair.first) << 32) | static_cast<size_t>(pair.second);
    }
};

void search(std::unordered_set<std::pair<int32_t, uint32_t>, pair_hash> & membership_result,
            std::vector<size_t> const & kmers,
            search_data const & data,
            search_config const & config,
            int64_t const ibf_idx)
{
    size_t const kmer_lemma = (kmers.size() > config.errors * config.k)
                              ? kmers.size() - config.errors * config.k
                              : 0;

    auto counting_agent = data.hibf[ibf_idx].template counting_agent<uint16_t>();

    auto && result = counting_agent.count_hashes(kmers);

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
