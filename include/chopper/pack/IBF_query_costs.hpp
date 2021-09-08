#pragma once

#include <cassert>
#include <stdexcept>

#include <robin_hood.h>

struct IBF_query_costs 
{
    static double get_exact(size_t const t_max)
    {
        // (19,19) mean c_{T_max} (see paper)
        static const robin_hood::unordered_map<uint32_t, double> costs =
        {
            {64, 1.0}, {128, 1.1}, {256, 1.32}, {512, 1.61},
            {1024, 2.69}, {2048, 4.48}, {4096, 7.53}, {8192, 13.65},
            {16384, 23.86}, {32768, 45.66}, {65536, 90.42}
        };

        if (costs.contains(t_max))
        {
            return costs.at(t_max);
        }
        else
        {
            throw std::invalid_argument("No exact data available for this t_max.");
        }
    }

    static double get_interpolated(size_t const t_max) 
    {
        if (t_max <= 64) return 1.0;

        // (19,19) mean c_{T_max} (see paper)
        static const robin_hood::unordered_map<uint32_t, double> costs =
        {
            {64, 1.0}, {128, 1.1}, {256, 1.32}, {512, 1.61},
            {1024, 2.69}, {2048, 4.48}, {4096, 7.53}, {8192, 13.65},
            {16384, 23.86}, {32768, 45.66}, {65536, 90.42}
        };
        
        if (costs.contains(t_max))
        {
            return costs.at(t_max);
        }

        for (size_t x = 128; x <= 65536; x *= 2)
        {
            if (t_max < x)
            {
                double const y0 = costs.at(x/2);
                double const y1 = costs.at(x);

                return y0 + (y1 - y0) * (t_max - x/2) / (x/2);
            }
        }

        throw std::invalid_argument("No data availabe for a t_max this large.");
    }
};