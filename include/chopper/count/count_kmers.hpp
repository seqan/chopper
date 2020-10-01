#include <future>
#include <iostream>
#include <thread>
#include <unordered_map>

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/range/views/to.hpp>

std::mutex mu;

void print_safely(std::pair<std::string, std::vector<std::string>> const & cluster,
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

void count_kmers(std::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                 count_config const & config)
{
    size_t const counting_threads = (config.num_threads <= 1) ? 1 : config.num_threads - 1;

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k}, seqan3::window_size{config.w});

    auto read_files = std::views::transform([] (auto const & cluster)
    {
        using result_t = std::vector<std::vector<seqan3::dna4>>;
        using pair_t = std::pair<std::pair<std::string, std::vector<std::string>>, result_t>;

        result_t result;
        for (auto const & filename : cluster.second)
        {
            sequence_file_type fin{filename};
            for (auto & [seq] : fin)
                result.push_back(seq);
        }

        return pair_t{cluster, result};
    });

    auto cluster_view = filename_clusters
                      | read_files
                      | seqan3::views::async_input_buffer(counting_threads);

    auto worker = [&cluster_view, &compute_minimiser] ()
    {
        for (auto && [cluster, sequence_vector] : cluster_view)
        {
            std::set<uint64_t> result;
            for (auto && seq : sequence_vector)
                for (auto hash : seq | compute_minimiser)
                    result.insert(hash);
            print_safely(cluster, result);
        }
    };

    // launch threads with worker
    std::vector<decltype(std::async(std::launch::async, worker))> handles;

    for (size_t i = 0; i < counting_threads; ++i)
        handles.emplace_back(std::async(std::launch::async, worker));
}
