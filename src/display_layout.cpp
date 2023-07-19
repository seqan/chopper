#include <iostream>
#include <set>

#include <robin_hood.h>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/stats/read_layout_file.hpp> // for read_layout_file

#include <hibf/detail/sketch/hyperloglog.hpp>

struct dna4_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<dna4_traits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

struct stats_configuration
{
    std::filesystem::path layout_file;
};

void keep_duplicates(std::vector<uint64_t> & shared, std::vector<uint64_t> const & current)
{
    auto shared_it = shared.begin();
    auto dup_it = shared.begin();
    auto current_it = current.begin();

    while (current_it != current.end() && shared_it != shared.end())
    {
        if (*current_it > *shared_it)
        {
            ++shared_it;
        }
        else if (*current_it == *shared_it)
        {
            std::swap(*dup_it, *shared_it);
            ++dup_it;
            ++shared_it;
            ++current_it;
        }
        else
        {
            ++current_it;
        }
    }

    shared.resize(dup_it - shared.begin());
}

int execute(stats_configuration & config)
{
    std::vector<std::vector<std::string>> filenames;
    chopper::configuration layout_config;
    hibf::layout::layout hibf_layout = chopper::stats::read_layout_file(layout_config, filenames, config.layout_file);

    std::sort(hibf_layout.user_bins.begin(), hibf_layout.user_bins.end());

    std::vector<hibf::sketch::hyperloglog> top_level_stats(layout_config.hibf_config.tmax,
                                                           layout_config.hibf_config.sketch_bits);
    std::vector<size_t> top_level_ub_counts(layout_config.hibf_config.tmax, 0);
    std::vector<std::vector<uint64_t>> top_level_shared_kmers(layout_config.hibf_config.tmax);
    std::vector<char> top_level_shared_kmers_state(layout_config.hibf_config.tmax, 'n');

    std::vector<uint64_t> current;

    bool const input_are_precomputed_files = filenames[0][0].ends_with(".minimiser");

    for (auto const & user_bin : hibf_layout.user_bins)
    {
        current.clear();

        size_t const idx =
            (user_bin.previous_TB_indices.size() == 0) ? user_bin.storage_TB_id : user_bin.previous_TB_indices[0];

        for (auto const & filename : filenames[user_bin.idx])
        {
            ++top_level_ub_counts[idx];

            if (input_are_precomputed_files)
            {
                uint64_t hash{};
                char * const hash_data{reinterpret_cast<char *>(&hash)};
                std::streamsize const hash_bytes{sizeof(hash)};

                std::ifstream infile{filename, std::ios::binary};

                while (infile.read(hash_data, hash_bytes))
                {
                    current.push_back(hash);
                    top_level_stats[idx].add(hash_data, hash_bytes);
                }
            }
            else
            {
                sequence_file_type fin{filename};

                for (auto && [seq] : fin)
                {
                    for (uint64_t hash_value : seq | seqan3::views::kmer_hash(seqan3::ungapped{layout_config.k}))
                    {
                        current.push_back(hash_value);
                        top_level_stats[idx].add(reinterpret_cast<char *>(&hash_value), sizeof(hash_value));
                    }
                }
            }

            std::sort(current.begin(), current.end());

            if (top_level_shared_kmers_state[idx] == 'n')
            {
                top_level_shared_kmers[idx] = current;
                top_level_shared_kmers_state[idx] = 'a';
            }
            else
            {
                keep_duplicates(top_level_shared_kmers[idx], current);
            }
        }
    }

    std::vector<chopper::layout::hibf_statistics::bin_kind> top_level_bin_kinds(
        layout_config.hibf_config.tmax,
        chopper::layout::hibf_statistics::bin_kind::split);

    for (auto const & max_bin : hibf_layout.max_bins)
    {
        if (max_bin.previous_TB_indices.size() == 1)
        {
            top_level_bin_kinds[max_bin.previous_TB_indices[0]] = chopper::layout::hibf_statistics::bin_kind::merged;
        }
    }

    // write out stats file
    std::cout << "tb_index\t"
              << "size\t"
              << "shared_size\t"
              << "ub_count\t"
              << "kind" << '\n';

    for (size_t i = 0; i < layout_config.hibf_config.tmax; ++i)
    {
        std::cout << i << '\t' << top_level_stats[i].estimate() << '\t' << top_level_shared_kmers[i].size() << '\t'
                  << top_level_ub_counts[i] << '\t'
                  << ((top_level_bin_kinds[i] == chopper::layout::hibf_statistics::bin_kind::merged) ? "merged"
                                                                                                     : "split")
                  << '\n';
    }

    return 0;
}

inline void set_up_parser(sharg::parser & parser, stats_configuration & config)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Compute an HIBF layout figure file";

    parser.info.description.emplace_back("Computes an table to siplay the layout.");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(
        config.layout_file,
        sharg::config{.short_id = '\0',
                      .long_id = "layout-file",
                      .description = "The input must be a layout file computed via chopper layout or raptor layout. ",
                      .required = true});
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"layout_stats", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    stats_configuration config;
    set_up_parser(parser, config);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    execute(config);
}
