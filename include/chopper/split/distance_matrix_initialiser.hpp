#pragma once

#include <chopper/split/minimizer.hpp>

struct distance_matrix_initialiser
{
    static constexpr size_t sketch_size{100};

    template <typename time_point>
    static std::string secs(time_point start, time_point end)
    {
        return "(" + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()) + "s)";
    }

    auto mash_distance(split_data & data, batch_config const & config)
    {
        map_distance_matrix distance_matrix{num_seq{length(data.sequences)},
                                            dummy_value{1.0},
                                            upper_distance_threshold{0.9}};

        // -----------------------------------------------------------------------------
        //                              SORT
        // -----------------------------------------------------------------------------
        auto start = std::chrono::steady_clock::now();

        // store the length s' < s of each sketch (some sequences are too small for full sketches)
        std::vector<size_t> sketch_lengths(length(data.sequences));

        // sort minimizer by value to get sketches
        for (size_t i = 0; i < length(data.sequences); ++i)
        {
            std::sort(seqan::begin(data.sequences[i]), seqan::end(data.sequences[i]));
            auto new_end_it = std::unique(seqan::begin(data.sequences[i]), seqan::end(data.sequences[i]));
            auto new_length = std::distance(seqan::begin(data.sequences[i]), new_end_it);
            sketch_lengths[i] = std::min<size_t>(sketch_size, new_length);
        }

        auto end = std::chrono::steady_clock::now();
        if (config.verbose)
            seqan3::debug_stream << ">>> Sorting individual minimizer sequences done. (" << secs(start, end) << std::endl;

        // -----------------------------------------------------------------------------
        // fill distance matrix
        // -----------------------------------------------------------------------------
        start = std::chrono::steady_clock::now();

        for (size_t i = 0; i < seqan::length(data.sequences); ++i)
        {
            for (size_t j = i + 1; j < seqan::length(data.sequences); ++j)
            {
                size_t common_hashes{0u};
                size_t unique_processed_hashes{0u};
                size_t pos_seq_i{0u};
                size_t pos_seq_j{0u};

                while (unique_processed_hashes < sketch_size &&
                       pos_seq_i < sketch_lengths[i] &&
                       pos_seq_j < sketch_lengths[j])
                {
                    minimizer const m_i = data.sequences[i][pos_seq_i];
                    minimizer const m_j = data.sequences[j][pos_seq_j];

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
                double const sim_score = (double)common_hashes / (double)unique_processed_hashes;
                distance_matrix.set_distance_value(i, j, similarity_score{sim_score}); // set to 1-score for distance
            }

            distance_matrix.set_distance_value(i, i, similarity_score{1.0}); // same sequence => similarity = 1
        }

        // std::cout << "distance matrices: ";
        // for (size_t i = 0; i < seqan::length(distance_matrix); ++i)
        //     std::cout << distance_matrix[i] << ',';
        // std::cout << std::endl;

        end = std::chrono::steady_clock::now();
        if (config.verbose)
            seqan3::debug_stream << ">>> Done filling distance matrix (" << secs(start, end) << std::endl;
        // -----------------------------------------------------------------------------

        // sort minimizer back
        for (size_t i = 0; i < length(data.sequences); ++i)
        {
            auto compare = [](minimizer const & m1, minimizer const & m2) { return m1.position < m2.position; };
            std::sort(seqan::begin(data.sequences[i]), seqan::end(data.sequences[i]), compare);
        }

        return distance_matrix;
    }
};
