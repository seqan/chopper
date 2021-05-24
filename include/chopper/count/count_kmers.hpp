#include <future>
#include <iostream>
#include <utility>
#include <thread>

#include <robin_hood.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/to.hpp>

std::mutex mu;

inline void print_safely(std::pair<std::string, std::vector<std::string>> const & cluster,
                         std::set<uint64_t> const & result)
{
    std::lock_guard<std::mutex> lock(mu);  // Acquire the mutex
    assert(cluster.second.size() >= 1);

    std::cout << cluster.second[0]; // print first filename
    for (size_t i = 1; i < cluster.second.size(); ++i)
        std::cout << ";" << cluster.second[i];
    std::cout << '\t' << result.size() << '\t' << cluster.first << std::endl;
}// lock_guard object is destroyed and mutex mu is released

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

template <typename cluster_view_type, typename compute_view_type>
void compute_hashes(cluster_view_type && cluster_view, compute_view_type && compute_fn)
{
    for (auto && [cluster, sequence_vector] : cluster_view)
    {
        std::set<uint64_t> result;
        for (auto && seq : sequence_vector)
            for (auto hash : seq | compute_fn)
                result.insert(hash);
        print_safely(cluster, result);
    }
}

inline void count_kmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                        count_config const & config)
{
    size_t const counting_threads = (config.num_threads <= 1) ? 1 : config.num_threads - 1;

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k}, seqan3::window_size{config.w});
    auto compute_kmers = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    auto read_files = std::views::transform([] (auto const & cluster)
    {
        using result_t = std::vector<std::vector<seqan3::dna4>>;
        using inner_pair_t = std::pair<std::string, std::vector<std::string>>;
        using pair_t = std::pair<inner_pair_t, result_t>;

        result_t result;
        for (auto const & filename : cluster.second)
        {
            sequence_file_type fin{filename};
            for (auto & [seq] : fin)
                result.push_back(seq);
        }

        return pair_t{inner_pair_t{cluster.first, cluster.second}, result};
    });

    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    auto filename = cluster_vector | seqan3::views::async_input_buffer(counting_threads * 5);

    auto worker = [&] ()
    {
        auto cluster_view = filename | read_files;

        if (config.disable_minimizers)
            compute_hashes(cluster_view, compute_kmers);
        else
            compute_hashes(cluster_view, compute_minimiser);
    };

    // launch threads with worker
    std::vector<decltype(std::async(std::launch::async, worker))> handles;

    for (size_t i = 0; i < counting_threads; ++i)
        handles.emplace_back(std::async(std::launch::async, worker));
}
