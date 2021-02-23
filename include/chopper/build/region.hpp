#pragma once

struct region
{
    size_t bin_index{};
    size_t begin{};
    size_t end{};

    bool operator==(region const & r) const
    {
        return std::tie(bin_index, begin, end) == std::tie(r.bin_index, r.begin, r.end);
    }

    bool operator!=(region const & r) const
    {
        return std::tie(bin_index, begin, end) != std::tie(r.bin_index, r.begin, r.end);
    }
};

inline std::ostream & operator<<(std::ostream & s, region const & r)
{
    s << "<" << r.bin_index << "," << r.begin << "," << r.end << ">";

    return s;
}
