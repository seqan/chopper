#pragma once

#include <thread>

#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::sketch
{

template <typename input_type>
    requires std::ranges::random_access_range<input_type> && /* input is a range of range of uint64_t values */
             std::ranges::input_range<std::ranges::range_value_t<input_type>>
          && std::same_as<uint64_t, std::ranges::range_value_t<std::ranges::range_value_t<input_type>>>
void compute_sketches_from_hashes(input_type && input,
                                  uint8_t const sketch_bits,
                                  uint8_t const number_of_threads,
                                  std::vector<chopper::sketch::hyperloglog> & sketches)
{
    size_t const input_size = std::ranges::size(input);
    sketches.resize(input_size);

#pragma omp parallel for schedule(static) num_threads(number_of_threads)
    for (size_t i = 0; i < input_size; ++i)
    {
        chopper::sketch::hyperloglog sketch(sketch_bits);

        for (auto && hash_value : input[i])
            sketch.add(reinterpret_cast<char *>(&hash_value), sizeof(hash_value));

#pragma omp critical
        sketches[i] = sketch;
    }
}

} // namespace chopper::sketch
