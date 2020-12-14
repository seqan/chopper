#pragma once

#include <chopper/split/map_distance_matrix.hpp>
#include <chopper/split/minimizer.hpp>
#include <chopper/split/split_config.hpp>
#include <chopper/split/split_data.hpp>

template <typename matrix_type = map_distance_matrix>
struct distance_matrix_initialiser
{
    static constexpr size_t sketch_size{100};

    size_t const number_of_sequences{};

    distance_matrix_initialiser(size_t const number_of_sequences_) :
        number_of_sequences{number_of_sequences_}
    {}

    template <typename time_point>
    static std::string secs(time_point start, time_point end)
    {
        return "(" + std::to_string(std::chrono::duration_cast<std::chrono::seconds>(end - start).count()) + "s)";
    }

    matrix_type initialise_matrix()
    {
        map_distance_matrix matrix{num_seq{number_of_sequences}, dummy_value{1.0}, upper_distance_threshold{0.9}};
        return matrix;
    }

    void set_distance_value(map_distance_matrix & matrix, size_t const i, size_t const j, double const dist)
    {
        matrix.set_distance_value(i, j, dist);
    }

    void set_distance_value(matrix_type & matrix, size_t const i, size_t const j, similarity_score const sim_score)
    {
        set_distance_value(matrix, i, j, 1.0 - sim_score.score);
    }

    void set_distance_value(matrix_type & matrix, size_t const i, size_t const j, distance_score const dis_score)
    {
        set_distance_value(matrix, i, j, dis_score.score);
    }

    matrix_type mash_distance(split_data & data, batch_config const & config)
    {
        assert(seqan::length(data.sequences) == number_of_sequences);
        matrix_type distance_matrix = initialise_matrix();

        // -----------------------------------------------------------------------------
        //                              SORT
        // -----------------------------------------------------------------------------
        auto start = std::chrono::steady_clock::now();

        // store the length s' < s of each sketch (some sequences are too small for full sketches)
        std::vector<size_t> sketch_lengths(number_of_sequences);

        // sort minimizer by value to get sketches
        for (size_t i = 0; i < number_of_sequences; ++i)
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

        for (size_t i = 0; i < number_of_sequences; ++i)
        {
            for (size_t j = i + 1; j < number_of_sequences; ++j)
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
                set_distance_value(distance_matrix, i, j, similarity_score{sim_score}); // set to 1-score for distance
            }

            set_distance_value(distance_matrix, i, i, similarity_score{1.0}); // same sequence => similarity = 1
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
        for (size_t i = 0; i < number_of_sequences; ++i)
        {
            auto compare = [](minimizer const & m1, minimizer const & m2) { return m1.position < m2.position; };
            std::sort(seqan::begin(data.sequences[i]), seqan::end(data.sequences[i]), compare);
        }

        return distance_matrix;
    }
};
