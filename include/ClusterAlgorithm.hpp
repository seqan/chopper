#pragma once

#include <seqan/graph_msa.h>
#include <seqan/seq_io.h>

#include "minimizer.hpp"
#include "segment_generation_config.hpp"

// globals
constexpr uint8_t kmer_size{25};
constexpr uint16_t window_size{100};

template <typename TNameSet>
bool
_loadSequences(seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
               TNameSet& fastaIDs,
               seqan::String<size_t> & original_sequence_lengths,
               const char *fileName)
{
    seqan::SeqFileIn inFile;
    if (!open(inFile, fileName))
    {
        std::cerr << "Could not open " << fileName << " for reading!" << std::endl;
        return false;
    }

    std::cerr << ">>> Processing file " << fileName << std::endl;
    seqan::StringSet<seqan::String<seqan::Dna>, seqan::Owner<>> sequences;
    {
        seqan::StringSet<seqan::String<seqan::Iupac>, seqan::Owner<>> iupac_sequences;
        seqan::readRecords(fastaIDs, iupac_sequences, inFile);
        seqan::resize(sequences, seqan::length(iupac_sequences));

        for (size_t idx = 0; idx < seqan::length(iupac_sequences); ++idx)
        {
            for (size_t jdx = 0; jdx < seqan::length(iupac_sequences[idx]); ++jdx)
            {
                seqan::appendValue(sequences[idx], iupac_sequences[idx][jdx]);
            }
            seqan::appendValue(original_sequence_lengths, seqan::length(iupac_sequences[idx]));
        }
    }

    // compute minimizers per sequence and store the corresponding chain in minimizer_sequences
    auto from = seqan::length(minimizer_sequences);
    resize(minimizer_sequences, seqan::length(minimizer_sequences) + seqan::length(sequences));
    Minimizer mini;
    mini.resize(kmer_size, window_size);
    for (size_t idx = from; idx < seqan::length(minimizer_sequences); ++idx)
    {
        minimizer_sequences[idx] = mini.getMinimizer(sequences[idx - from]);
        // seqan3::debug_stream << seqan::length(minimizer_sequences[idx]) << std::endl;
    }

    return (seqan::length(fastaIDs) > 0u);
}

struct ClusterAlgorithm
{
    static constexpr double FPR{1.0E-6};
    static constexpr size_t sketch_size{100};
    static constexpr uint8_t hash_funs{2}; // nuber of hash functions
    static constexpr float similarity_threshold{0.8f};
    static constexpr float containment_index_threshold{0.2f};

    using fitting_bin_info_t = std::vector<std::pair<size_t, seqan3::bin_index>>;
    using fitting_bin_t = std::vector<std::vector<fitting_bin_info_t>>;

    template <typename time_point>
    std::string secs(time_point start, time_point end)
    {
        return "(" + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()) + "s)";
    }


    size_t map_length_to_ibf_pos(size_t const length, std::vector<size_t> const & ranges)
    {
        assert(std::lower_bound(ranges.begin(), ranges.end(), length) != ranges.end());
        return std::distance(ranges.begin(), std::lower_bound(ranges.begin(), ranges.end(), length));
    }

    template <typename begin_iterator_type, typename end_iterator_type>
    auto get_kmer_counts(begin_iterator_type const sketch_begin,
                         end_iterator_type const sketch_end,
                         seqan3::interleaved_bloom_filter<> const & ibf)
    {
        auto agent = ibf.membership_agent();
        assert(ibf.bin_count() % 64 == 0);
        std::vector<size_t> counts(ibf.bin_count() / 64);
        for (auto it = sketch_begin; it != sketch_end; ++it)
            for (auto && [idx, contained] : seqan3::views::zip(std::views::iota(0), agent.bulk_contains(*it)))
                counts[idx / 64] += contained; // contained = 0/1
        return counts;
    }

    template <typename begin_iterator_type, typename end_iterator_type>
    seqan3::bin_index there_exists_a_bin_that_fits(begin_iterator_type const sketch_begin,
                                                   end_iterator_type const sketch_end,
                                                   seqan3::interleaved_bloom_filter<> const & ibf)
    {
        size_t const sz = std::distance(sketch_begin, sketch_end);
        auto const counts = get_kmer_counts(sketch_begin, sketch_end, ibf);

        auto max_it = std::max_element(counts.begin(), counts.end());

        if ((double)*max_it / (double)sz > similarity_threshold)
            return seqan3::bin_index{(size_t)std::distance(counts.begin(), max_it) * 64};
        else
            return seqan3::bin_index{std::numeric_limits<size_t>::max()};
    }

    template <typename begin_iterator_type, typename end_iterator_type>
    void fill_64_bins_with(seqan3::interleaved_bloom_filter<> & ibf,
                           size_t const start_bin,
                           begin_iterator_type const begin_it,
                           end_iterator_type const end_it)
    {
        // add all kmers of seq to that bin
        // chunk seq in 64 bins
        size_t const sz = std::distance(begin_it, end_it);
        auto const chunk_size{(sz + 64)/64};

        for (size_t chunk_idx = 0; chunk_idx * chunk_size < sz; ++chunk_idx)
        {
            seqan3::bin_index const bidx{start_bin + chunk_idx};
            auto const start = begin_it + chunk_size * chunk_idx;
            auto const stop = begin_it + std::min<size_t>(chunk_size * (chunk_idx + 1), sz);

            for (auto it = start; it != stop; ++it)
                ibf.emplace(*it, bidx);
        }
    }

    void fill_ranges_and_bin_sizes(std::vector<size_t> & ranges,
                                   std::vector<size_t> & bin_sizes,
                                   seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences)
    {
        auto length_comparator = [] (auto const & seq1, auto const & seq2) { return length(seq1) < length(seq2); };
        auto const [minit, maxit] = std::minmax_element(seqan::begin(minimizer_sequences),
                                                        seqan::end(minimizer_sequences),
                                                        length_comparator);

        auto const min = std::max<size_t>(seqan::length(*minit), 129u); // Too small minimizer sequences have no point for an IBF
        auto const max = std::max<size_t>(seqan::length(*maxit), min);

        size_t const num_ibfs = std::floor((std::log(min) - std::log(max))/std::log(2.0/3.0)) + 1;

    seqan3::debug_stream << "    min: " << min << " max: " << max << " num_ibfs: " << num_ibfs << std::endl;

        ranges.resize(num_ibfs);
        bin_sizes.resize(num_ibfs);

        ranges[num_ibfs - 1] = max + 1; // plus 1 s.t. that std::lower_bound returns the last element for the max
        for (int i = num_ibfs - 2; i >= 0; --i)
            ranges[i] = ranges[i + 1] * 2 / 3;

        for (size_t ibf_idx = 0; ibf_idx < num_ibfs; ++ibf_idx)
        {
            auto const n = ranges[ibf_idx]; // take the maximum as bloom filter fill
            auto const m = (size_t)std::ceil((-(double)(hash_funs * ((n + 64)/64)) / // /64 because one sequence always gets 64 bins
                                             (std::log(1 - std::pow(10.0, std::log10(FPR) / hash_funs)))));
            bin_sizes[ibf_idx] = m;
        }
    }

    template <typename TNameSet>
    void sort_by_last_parameters_length(TNameSet & names,
                                        seqan::String<size_t> & original_sequence_lengths,
                                        seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & sequences)
    {
        // generate permutation of indices sorted in descinding order by the sequence lengths
        auto permutation = std::views::iota(0u, length(sequences)) | seqan3::views::to<std::vector>;
        assert(permutation.size() == length(sequences));
        auto index_compare = [&sequences] (auto const i1, auto const i2)
                             { return length(sequences[i2]) < length(sequences[i1]); };
        std::sort(permutation.begin(), permutation.end(), index_compare);

        // apply permutation on names and sequences
        for (size_t i = 0; i < permutation.size(); i++)
        {
            auto current = i;
            while (i != permutation[current])
            {
                auto next = permutation[current];
                std::swap(names[current], names[next]);
                std::swap(sequences[current], sequences[next]);
                std::swap(original_sequence_lengths[current], original_sequence_lengths[next]);
                permutation[current] = current;
                current = next;
            }
            permutation[current] = current;
        }
    }

    template <typename begin_iterator_type, typename end_iterator_type>
    void add_fitting_bins(fitting_bin_t & seq_fitting_bins,
                          size_t const ibf_idx,
                          size_t const ibf_bin_pos,
                          begin_iterator_type const sketch_begin,
                          end_iterator_type const sketch_end,
                          std::vector<std::optional<seqan3::interleaved_bloom_filter<>>> const & ibfs)
    {
        seq_fitting_bins[ibf_idx].push_back(fitting_bin_info_t{}); // init
        size_t const sz = std::distance(sketch_begin, sketch_end);

        for (size_t idx = ibf_idx + 1; idx < ibfs.size(); ++idx)
        {
            if (ibfs[idx].has_value())
            {
                auto const counts = get_kmer_counts(sketch_begin, sketch_end, ibfs[idx].value());

                for (size_t i = 0; i < counts.size(); ++i)
                    if ((double)counts[i] / (double)sz > containment_index_threshold)
                        seq_fitting_bins[ibf_idx][ibf_bin_pos].emplace_back(idx, std::move(seqan3::bin_index{i}));
            }

        }
    }

    template<typename TSize2, typename TSpec2>
    inline void add_all_to_all(std::vector<size_t> const & seq_ids,
                               seqan::String<TSize2, TSpec2>& pList)
    {
        // reserve(pList, length(pList) + seq_ids_small.size() * seq_ids_large.size());

        for (size_t i = 0; i < seq_ids.size(); ++i)
        {
            for (size_t j = i + 1; j < seq_ids.size(); ++j)
            {
                seqan::appendValue(pList, seq_ids[i]);
                seqan::appendValue(pList, seq_ids[j]);
            }
        }
    }

    template<typename TSize2, typename TSpec2>
    inline void add_all_to_first(std::vector<size_t> const & seq_ids,
                                 seqan::String<TSize2, TSpec2>& pList)
    {
        for (size_t i = 1; i < seq_ids.size(); ++i)
        {
            seqan::appendValue(pList, seq_ids[0]);
            seqan::appendValue(pList, seq_ids[i]);
        }
    }

    template<typename TSize2, typename TSpec2>
    inline void add_many_to_many(std::vector<size_t> const & seq_ids_small,
                                 std::vector<size_t> const & seq_ids_large,
                                 seqan::String<TSize2, TSpec2>& pList)
    {
        // reserve(pList, length(pList) + seq_ids_small.size() * seq_ids_large.size());

        for (size_t i = 0; i < seq_ids_small.size(); ++i)
        {
            for (size_t j = 0; j < seq_ids_large.size(); ++j)
            {
                seqan::appendValue(pList, seq_ids_small[i]);
                seqan::appendValue(pList, seq_ids_large[j]);
            }
        }
    }

    void set_distance_value(seqan::String<double> & distanceMatrix, size_t i, size_t j, double similarity_score)
    {
        size_t const nseq = std::sqrt(length(distanceMatrix));

        // since we are not sure which diagonal is used but they are symmetric, fill both diagonals
        if (i < j)
            distanceMatrix[i * nseq + j] = 1.0 - similarity_score; // distance = 1 - similarity
        else
            distanceMatrix[j * nseq + i] = 1.0 - similarity_score; // distance = 1 - similarity
    }

    void set_distance_value(map_distance_matrix & distanceMatrix, size_t i, size_t j, double similarity_score)
    {
        size_t const nseq = std::sqrt(seqan::length(distanceMatrix));

        // always set whole upper and lower diagonal because I don't know what is needed
        distanceMatrix.emplace(i * nseq + j, 1.0 - similarity_score); // distance = 1 - similarity
        distanceMatrix.emplace(j * nseq + i, 1.0 - similarity_score); // distance = 1 - similarity
    }

    void set_all_to_all_to(seqan::String<double> & distanceMatrix, std::vector<size_t> const & ids, double score)
    {
        for (size_t i = 0; i < ids.size(); ++i)
        {
            for (size_t j = i + 1; j < ids.size(); ++j)
            {
                set_distance_value(distanceMatrix, ids[i], ids[j], score);
            }
        }
    }

    void set_one_to_many_to(seqan::String<double> & distanceMatrix, size_t i, std::vector<size_t> const & ids, double score)
    {
        for (size_t j = 0; j < ids.size(); ++j)
        {
            set_distance_value(distanceMatrix, i, ids[j], score);
        }
    }

    void fill_distance_matrix(seqan::String<double> & distanceMatrix,
                              seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
                              std::vector<std::vector<std::vector<size_t>>> const & ibfs_seq_indices,
                              std::vector<std::optional<seqan3::interleaved_bloom_filter<>>> const & ibfs)
    {
        auto const nseq = length(minimizer_sequences);
        resize(distanceMatrix, nseq * nseq, 0.0);

        // for (size_t i = 0; i < nseq; ++i)
        //     set_distance_value(distanceMatrix, i, i, 0); // set diagonal to 0 distance

        // for each ibf
        for (size_t ibf_idx = 0; ibf_idx < ibfs_seq_indices.size(); ++ibf_idx)
        {
            // for each bin in each ibf
            for (size_t bin_idx = 0; bin_idx < ibfs_seq_indices[ibf_idx].size(); ++bin_idx)
            {
                // the ids of this bin all have a similarity of > 0.9 to the first entry in ids
                auto const & ids = ibfs_seq_indices[ibf_idx][bin_idx];
                set_all_to_all_to(distanceMatrix, ids, 0.9);

                // add the similarity between bins of the same ibf
                for (size_t i = 0; i < ids.size(); ++i)
                {
                    auto const current_seq_idx = ids[i];
                    auto const seq_begin = seqan::begin(minimizer_sequences[current_seq_idx]);
                    auto const new_end_it = std::unique(seq_begin, seqan::end(minimizer_sequences[current_seq_idx]));
                    auto const sketch_end = std::min(seq_begin + sketch_size, new_end_it);
                    auto const sz = std::distance(seq_begin, sketch_end);

                    auto const counts = get_kmer_counts(seq_begin, sketch_end, ibfs[ibf_idx].value());

                    assert(counts.size() == ibfs_seq_indices[ibf_idx].size());

                    for (size_t other_bin_idx = bin_idx; other_bin_idx < counts.size(); ++other_bin_idx)
                    {
                        auto const & other_ids = ibfs_seq_indices[ibf_idx][other_bin_idx];
                        double const score = (double)counts[other_bin_idx] / (double)sz;

                        set_one_to_many_to(distanceMatrix, current_seq_idx, other_ids, score);
                    }
                }

                // add the similarity between bins of a greater ibf
                for (size_t other_ibf_idx = ibf_idx + 1; other_ibf_idx < ibfs_seq_indices.size(); ++other_ibf_idx)
                {
                    if (ibfs[other_ibf_idx].has_value())
                    {
                        for (size_t i = 0; i < ids.size(); ++i)
                        {
                            auto const current_seq_idx = ids[i];
                            auto const seq_begin = seqan::begin(minimizer_sequences[current_seq_idx]);
                            auto const new_end_it = std::unique(seq_begin, seqan::end(minimizer_sequences[current_seq_idx]));
                            auto const sketch_end = std::min(seq_begin + sketch_size, new_end_it);
                            auto const sz = std::distance(seq_begin, sketch_end);

                            auto const counts = get_kmer_counts(seq_begin, sketch_end, ibfs[other_ibf_idx].value());

                            assert(counts.size() == ibfs_seq_indices[other_ibf_idx].size());

                            for (size_t other_bin_idx = 0; other_bin_idx < counts.size(); ++other_bin_idx)
                            {
                                auto const & other_ids = ibfs_seq_indices[other_ibf_idx][other_bin_idx];
                                double const size_ratio = seqan::length(minimizer_sequences[current_seq_idx]) /
                                                          seqan::length(minimizer_sequences[ibfs_seq_indices[other_ibf_idx][other_bin_idx][0]]);
                                double const score = ((double)counts[other_bin_idx] / (double)sz) * size_ratio; // scale similarity by size difference

                                set_one_to_many_to(distanceMatrix, current_seq_idx, other_ids, score);
                            }
                        }
                    }
                }
            }
        }
    }

    template <typename TNameSet, typename TSize>
    auto load_sequences_of_multiple_files(seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
                                          TNameSet & fastaIDs,
                                          seqan::String<size_t> & original_sequence_lengths,
                                          segment_generation_config<TSize> & config)
    {
        // -----------------------------------------------------------------------------
        //                              LOAD DATA
        // -----------------------------------------------------------------------------
        auto start = std::chrono::steady_clock::now();

            for (auto const & file_name : config.seqfiles)
                if (!_loadSequences(minimizer_sequences, fastaIDs, original_sequence_lengths, file_name.c_str()))
                    return false;

        auto end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Loading sequences and computing minimizers complete " << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------
        //                     FILL RANGES AND BIN_SIZES
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

            std::vector<size_t> ranges;
            std::vector<size_t> bin_sizes;
            fill_ranges_and_bin_sizes(ranges, bin_sizes, minimizer_sequences);

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Computing ranges done "  << secs(start, end)
                             << "    ranges: " << ranges << std::endl
                             << "    IBF bin sizes: " << bin_sizes << std::endl;
        // -----------------------------------------------------------------------------
        //                              SORT
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

            // sort minimizer sequences in descending order
            sort_by_last_parameters_length(fastaIDs, original_sequence_lengths, minimizer_sequences);

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Sorting Done. (" << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------
        //                              CLUSTER
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

            std::vector<std::optional<seqan3::interleaved_bloom_filter<>>> ibfs(ranges.size());
            std::vector<std::vector<std::vector<size_t>>> ibfs_seq_indices(ranges.size());

            fitting_bin_t seq_fitting_bins(ranges.size());

            for (size_t i = 0; i < length(minimizer_sequences); ++i)
            {
                std::sort(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]));

                auto new_end_it = std::unique(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]));
                auto l = std::distance(seqan::begin(minimizer_sequences[i]), new_end_it);

                auto ibf_idx = map_length_to_ibf_pos(l, ranges);
                auto const seq_begin = seqan::begin(minimizer_sequences[i]);
                auto const sketch_end = std::min(seq_begin + sketch_size, new_end_it);

        // seqan3::debug_stream << "    " << fastaIDs[i];
        // seqan3::debug_stream << "\t(" << std::distance(new_end_it, seqan::end(minimizer_sequences[i])) << " duplicates) ";
        // seqan3::debug_stream << "\tlen:" << l << "\tibf:" << ibf_idx;

                if (!ibfs[ibf_idx].has_value()) // ibf was not initialised yet
                {
                    ibfs[ibf_idx] = seqan3::interleaved_bloom_filter{seqan3::bin_count{64},
                                                                     seqan3::bin_size{bin_sizes[ibf_idx]},
                                                                     seqan3::hash_function_count{hash_funs}};
                    ibfs_seq_indices[ibf_idx].push_back(std::vector<size_t>{i});
                    fill_64_bins_with(ibfs[ibf_idx].value(), 0, seq_begin, new_end_it);

                    // check ibfs that contain larger sequences.
                    // because we have sorted our sequences in descending order, we do not need to check smaller ones
                    add_fitting_bins(seq_fitting_bins, ibf_idx, 0, seq_begin, sketch_end, ibfs);
        // seqan3::debug_stream << "\tbin_idx:" << 0;
                }
                else if (seqan3::bin_index bin = there_exists_a_bin_that_fits(seq_begin, sketch_end, ibfs[ibf_idx].value());
                         bin.get() != std::numeric_limits<size_t>::max())
                {
                    ibfs_seq_indices[ibf_idx][(bin.get()/64)].push_back(i);
        // seqan3::debug_stream << "\tbin_idx:" << (bin.get()/64);
                }
                else
                {
                    // add a new bin to the ibf of the respective sizes
                    size_t previous_ibf_size{ibfs[ibf_idx]->bin_count()};
                    ibfs[ibf_idx]->increase_bin_number_to(seqan3::bin_count{previous_ibf_size + 64});
                    ibfs_seq_indices[ibf_idx].push_back(std::vector<size_t>{i});
                    fill_64_bins_with(ibfs[ibf_idx].value(), previous_ibf_size, seq_begin, new_end_it);

                    // check ibfs that contain larger sequences.
                    // because we have sorted our sequences in descending order, we do not need to check smaller ones
                    add_fitting_bins(seq_fitting_bins, ibf_idx, (previous_ibf_size/64), seq_begin, sketch_end, ibfs);
        // seqan3::debug_stream << "\tbin_idx:" << (previous_ibf_size/60);
        // seqan3::debug_stream << "\tfitting bins [";
        // for (auto & p : seq_fitting_bins[ibf_idx][(previous_ibf_size/60)])
        //     seqan3::debug_stream  << "("<< p.first << "," << p.second.get() << "),";
        // seqan3::debug_stream << "]";
                }
        // seqan3::debug_stream << std::endl;
            }

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Done clustering (" << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------

        // -----------------------------------------------------------------------------
        // fill distance matrix
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

        fill_distance_matrix(config.distanceMatrix, minimizer_sequences, ibfs_seq_indices, ibfs);

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Done filling distance matrix (" << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------

        // sort minimizer back
        for (size_t i = 0; i < length(minimizer_sequences); ++i)
        {
            auto compare = [](minimizer const & m1, minimizer const & m2) { return m1.position < m2.position; };
            std::sort(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]), compare);
        }

        for (auto const & [per_ibf_ids, per_ibf_fitting_bins] : seqan3::views::zip(ibfs_seq_indices, seq_fitting_bins))
        {
            for (auto const & [ids, fitting_bins] : seqan3::views::zip(per_ibf_ids, per_ibf_fitting_bins))
            {
                // all to all pairwise global alignments within the same ibf
                add_all_to_all(ids, config.global_alignment_pairs);

                for (auto const & [ibf_idx, bin_idx] : fitting_bins)
                    add_many_to_many(ids, ibfs_seq_indices[ibf_idx][bin_idx.get()], config.semi_global_alignment_pairs);
            }

            auto joined_ids = seqan3::views::join(per_ibf_ids) | seqan3::views::to<std::vector>;
            add_all_to_all(joined_ids, config.local_alignment_pairs);
        }

        return true;
    }

    template <typename TNameSet, typename TSize>
    auto fill_mash_distance_matrix(seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
                                   TNameSet & fastaIDs,
                                   seqan::String<size_t> & original_sequence_lengths,
                                   segment_generation_config<TSize> & config)
    {
        // -----------------------------------------------------------------------------
        //                              LOAD DATA
        // -----------------------------------------------------------------------------
        auto start = std::chrono::steady_clock::now();

        for (auto const & file_name : config.seqfiles)
            if (!_loadSequences(minimizer_sequences, fastaIDs, original_sequence_lengths, file_name.c_str()))
                return false;

        auto end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Loading sequences and computing minimizers complete " << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------
        //                              SORT
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

        // store the length s' < s of each sketch (some sequences are too small for full sketches)
        std::vector<size_t> sketch_lengths(length(minimizer_sequences));

        // sort minimizer by value to get sketches
        for (size_t i = 0; i < length(minimizer_sequences); ++i)
        {
            std::sort(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]));
            auto new_end_it = std::unique(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]));
            auto new_length = std::distance(seqan::begin(minimizer_sequences[i]), new_end_it);
            sketch_lengths[i] = std::min<size_t>(sketch_size, new_length);
        }

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Sorting individual minimizer sequences done. (" << secs(start, end) << std::endl;

        // -----------------------------------------------------------------------------
        // fill distance matrix
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

        config.distanceMatrix.nseq = length(minimizer_sequences);

        for (size_t i = 0; i < seqan::length(minimizer_sequences); ++i)
        {
            for (size_t j = i + 1; j < seqan::length(minimizer_sequences); ++j)
            {
                size_t common_hashes{0u};
                size_t unique_processed_hashes{0u};
                size_t pos_seq_i{0u};
                size_t pos_seq_j{0u};

                while (unique_processed_hashes < sketch_size &&
                       pos_seq_i < sketch_lengths[i] &&
                       pos_seq_j < sketch_lengths[j])
                {
                    minimizer const m_i = minimizer_sequences[i][pos_seq_i];
                    minimizer const m_j = minimizer_sequences[j][pos_seq_j];

                    pos_seq_i     += m_i <= m_j;
                    pos_seq_j     += m_i >= m_j;
                    common_hashes += m_i == m_j;

                    ++unique_processed_hashes;
                }

                // to avoid an if clause here, we always do this small computation
                // check performance if "if (unique_processed_hashes < sketch_size)" is more efficient.
                size_t const remaining = sketch_size - unique_processed_hashes;
                size_t const available_i = sketch_lengths[i] - pos_seq_i;
                size_t const available_j = sketch_lengths[j] - pos_seq_j;
                assert(remaining == 0 || available_i == 0 || available_j == 0);
                unique_processed_hashes += std::min(remaining, available_i + available_j);

                // Jaquard Index is approximated with x/s' -> common_hashes / unique_processed_hashes
                double const similarity_score = (double)common_hashes / (double)unique_processed_hashes;
                if (similarity_score > 0.1)
                    set_distance_value(config.distanceMatrix, i, j, similarity_score); // set to 1-score for distance
            }

            set_distance_value(config.distanceMatrix, i, i, 1); // same sequence => similarity = 1
        }

        std::cout << "distance matrices: ";
        for (size_t i = 0; i < seqan::length(config.distanceMatrix); ++i)
            std::cout << config.distanceMatrix[i] << ',';
        std::cout << std::endl;

        end = std::chrono::steady_clock::now();
        seqan3::debug_stream << ">>> Done filling distance matrix (" << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------

        // sort minimizer back
        for (size_t i = 0; i < length(minimizer_sequences); ++i)
        {
            auto compare = [](minimizer const & m1, minimizer const & m2) { return m1.position < m2.position; };
            std::sort(seqan::begin(minimizer_sequences[i]), seqan::end(minimizer_sequences[i]), compare);
        }

        return true;
    }
};
