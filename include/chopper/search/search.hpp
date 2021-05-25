#pragma once

#include <chopper/search/pair_hash.hpp>
#include <chopper/search/search_config.hpp>
#include <chopper/search/search_data.hpp>

inline void clear_and_compute_kmers(std::vector<size_t> & kmers, seqan3::dna4_vector const & query, search_config const & config)
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

inline void write_result(std::vector<std::pair<int32_t, uint32_t>> & membership_result,
                         std::string const & id,
                         search_data const & data,
                         std::ostream & out_stream)
{
    out_stream << id << '\t';

    if (membership_result.empty())
    {
        out_stream << std::endl;
        return;
    }

    // otherwise the result output is not testable
    std::ranges::sort(membership_result, [&data] (auto const & pair1, auto const & pair2)
                                         {
                                             return data.user_bins.filename_index(pair1.first, pair1.second) <
                                                    data.user_bins.filename_index(pair2.first, pair2.second);
                                         });

    for (size_t i = 0; i < membership_result.size() - 1; ++i)
    {
        auto & [ibf_idx, bin_idx] = membership_result[i];
        assert(data.user_bins.filename_index(ibf_idx, bin_idx) > -1);
        out_stream << data.user_bins.filename_index(ibf_idx, bin_idx) << ',';
    }

    out_stream << data.user_bins.filename_index(membership_result.back().first, membership_result.back().second) << std::endl;
}

inline void search(std::vector<std::pair<int32_t, uint32_t>> & membership_result,
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
    assert(result.size() > 0);

    size_t sum{};

    for (size_t bin{}; bin < result.size() - 1; ++bin)
    {
        sum += result[bin];
        auto const current_filename_index = data.user_bins.filename_index(ibf_idx, bin);

        if (current_filename_index < 0) // merged bin
        {
            // if threshold, next level
            if (sum >= kmer_lemma)
                search(membership_result, kmers, data, config, data.hibf_bin_levels[ibf_idx][bin]);
            sum = 0;
        }
        else if (current_filename_index != data.user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
        {
            // if threshold, write
            if (sum >= kmer_lemma)
                membership_result.emplace_back(ibf_idx, bin);
            sum = 0;
        }
    }

    // check the last bin
    if (sum + result.back() >= kmer_lemma)
        if (auto bin =  result.size() - 1; data.user_bins.filename_index(ibf_idx, bin) < 0)
            search(membership_result, kmers, data, config, data.hibf_bin_levels[ibf_idx][bin]);
        else
            membership_result.emplace_back(ibf_idx, bin);
}
