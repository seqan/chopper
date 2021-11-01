#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <future>
#include <thread>

#include <robin_hood.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/to.hpp>

#include <chopper/count/configuration.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::count
{

inline void write_cluster_data(std::pair<std::string, std::vector<std::string>> const & cluster,
                               uint64_t const size,
                               std::ofstream & fout)
{
    assert(cluster.second.size() >= 1);

    fout << cluster.second[0]; // write first filename
    for (size_t i = 1; i < cluster.second.size(); ++i)
        fout << ";" << cluster.second[i];
    fout << '\t' << size << '\t' << cluster.first << std::endl;
}

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

template <typename seq_type, typename compute_view_type>
void compute_hashes(seq_type && seq,
                    compute_view_type && compute_fn,
                    configuration const & config,
                    robin_hood::unordered_node_set<uint64_t> & result,
                    chopper::sketch::hyperloglog & sketch)
{
    if (!config.exclusively_hlls)
    {
        auto hash_view = seq | compute_fn | std::views::common;
        result.insert(hash_view.begin(), hash_view.end());
    }

    if (config.exclusively_hlls || !config.hll_dir.empty())
        for (auto && hash : seq | compute_fn)
            sketch.add(reinterpret_cast<char*>(&hash), sizeof(hash));
}

inline void count_kmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                        configuration const & config)
{
    // output file
    std::ofstream fout{config.output_filename};

    // create the hll dir if it doesn't already exist
    if (!config.hll_dir.empty())
        std::filesystem::create_directory(config.hll_dir);

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k},
                                                           seqan3::window_size{config.w},
                                                           seqan3::seed{0x8F3F73B5CF1C9ADE >> (64u - 2u * config.k)});
    auto compute_kmers = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    #pragma omp parallel for schedule(static) num_threads(config.num_threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        // read files
        std::vector<std::vector<seqan3::dna4>> sequence_vector{};
        for (auto const & filename : cluster_vector[i].second)
            for (auto && [seq] : sequence_file_type{filename})
                sequence_vector.push_back(seq);

        robin_hood::unordered_node_set<uint64_t> result{};
        chopper::sketch::hyperloglog sketch(config.sketch_bits);

        if (config.disable_minimizers)
            for (auto && seq : sequence_vector)
                compute_hashes(seq, compute_kmers, config, result, sketch);
        else
            for (auto && seq : sequence_vector)
                compute_hashes(seq, compute_minimiser, config, result, sketch);

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t const size = config.exclusively_hlls ? static_cast<uint64_t>(sketch.estimate()) : result.size();

        #pragma omp critical
        write_cluster_data(cluster_vector[i], size, fout);

        if (!config.hll_dir.empty())
        {
            // For more than one file in the cluster, Felix doesn't know how to name the file
            // and what exactly is supposed to happen.
            if (cluster_vector[i].second.size() > 1)
                throw std::runtime_error("This mode is not implemented yet for multiple files grouped together.");

            // For one file in the cluster, the file stem is used with the .hll ending
            std::filesystem::path path = config.hll_dir / std::filesystem::path(cluster_vector[i].first).stem();
            path += ".hll";
            std::ofstream hll_fout(path, std::ios::binary);
            sketch.dump(hll_fout);
        }
    }
}

} // namespace chopper::count
