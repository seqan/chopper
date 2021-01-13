#pragma once

struct region
{
    size_t hidx{};
    int64_t lidx{};
    size_t begin{};
    size_t end{};

    bool operator==(region const & r) const
    {
        return std::tie(hidx, lidx, begin, end) == std::tie(r.hidx, r.lidx, r.begin, r.end);
    }

    bool operator!=(region const & r) const
    {
        return std::tie(hidx, lidx, begin, end) != std::tie(r.hidx, r.lidx, r.begin, r.end);
    }
};

std::ostream & operator<<(std::ostream & s, region const & r)
{
    s << "<" << r.hidx << "," << r.lidx << "," << r.begin << "," << r.end << ">";

    return s;
}
