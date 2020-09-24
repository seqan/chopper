#pragma once

#include "minimizer.hpp"
#include "segment_generation_config.hpp"

struct ClusterAlgorithm
{
    static constexpr size_t sketch_size{100};

    template <typename time_point>
    static std::string secs(time_point start, time_point end)
    {
        return "(" + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()) + "s)";
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

    template <typename TNameSet, typename TSize>
    auto fill_mash_distance_matrix(seqan::StringSet<seqan::String<minimizer>, seqan::Owner<>> & minimizer_sequences,
                                   TNameSet & fastaIDs,
                                   seqan::String<size_t> & original_sequence_lengths,
                                   segment_generation_config<TSize> & config)
    {
        // -----------------------------------------------------------------------------
        //                              SORT
        // -----------------------------------------------------------------------------
        auto start = std::chrono::steady_clock::now();

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

        auto end = std::chrono::steady_clock::now();
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

        // std::cout << "distance matrices: ";
        // for (size_t i = 0; i < seqan::length(config.distanceMatrix); ++i)
        //     std::cout << config.distanceMatrix[i] << ',';
        // std::cout << std::endl;

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
