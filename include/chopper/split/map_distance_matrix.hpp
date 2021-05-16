#pragma once

#include <robin_hood.h>

struct similarity_score
{
    double score;
};

struct distance_score
{
    double score;
};


struct num_seq
{
    size_t n;
};

struct dummy_value
{
    double v;
};

struct upper_distance_threshold
{
    double t;
};

struct map_distance_matrix : robin_hood::unordered_map<size_t, double>
{
    using base_t = robin_hood::unordered_map<size_t, double>;
    using value_type = double;
    using size_type = size_t;

    struct proxy
    {
        /*!\name Construction, destruction and assignment
         * \{
         */
        proxy() = default; //!\brief Default construction.
        proxy(proxy const & rhs) = default; //!\brief Copy construction.
        proxy(proxy && rhs) = default; //!\brief Move construction.
        proxy & operator=(proxy const & rhs) = default; //!\brief Copy assignment.
        proxy & operator=(proxy && rhs) = default; //!\brief Move assignment.
        ~proxy() = default; //!\brief Destruction.

        //!\brief Constructing from host and value
        proxy(map_distance_matrix & matrix, double val, size_t index_) :
            host{&matrix}, value{std::move(val)}, index{std::move(index_)}
        {}

        // assign from value type
        proxy & operator=(value_type const & val)
        {
            host->set_distance_value(index, val);
            value = host->get_distance_value(index);
            return *this;
        }

        proxy & operator=(value_type && val)
        {
            host->set_distance_value(index, val);
            value = host->get_distance_value(index);
            return *this;
        } //!\brief Move assignment.
        //!\}

        operator value_type() const
        {
            return value;
        }

    private:
        map_distance_matrix * host{nullptr};
        double value{};
        size_t index{};
    };

    map_distance_matrix(num_seq const number_of_sequences,
                        dummy_value const dummy,
                        upper_distance_threshold const threshold) :
        nseq{number_of_sequences.n},
        dummy{dummy.v},
        distance_threshold{threshold.t}
    {}

    size_type const nseq{};

    value_type const dummy{};

    value_type const distance_threshold{};


    proxy operator[](size_type index)
    {
        return proxy{*this, get_distance_value(index), index};
    }

    // i = first = row , j = second = column
    // 0 1 2 3 4     i * nseq + j
    // 5 6 7 8 9
    //
    void set_distance_value(size_t i, size_t j, similarity_score sim_score)
    {
        assert(nseq != 0);
        set_distance_value(i * nseq + j, 1.0 - sim_score.score);
        set_distance_value(j * nseq + i, 1.0 - sim_score.score);
    }

    void set_distance_value(size_t i, size_t j, distance_score dis_score)
    {
        assert(nseq != 0);
        set_distance_value(i * nseq + j, dis_score.score);
        set_distance_value(j * nseq + i, dis_score.score);
    }

    void set_distance_value(size_t i, size_t j, double score)
    {
        assert(nseq != 0);
        set_distance_value(i * nseq + j, score);
        set_distance_value(j * nseq + i, score);
    }

// private:
    void set_distance_value(size_t index, value_type distance_value)
    {
        assert(nseq != 0);

        auto it = this->find(index);

        if (it == this->end()) // index has not been inserted yet
        {
            if (distance_value <= distance_threshold)
                this->emplace(index, distance_value);
        }
        else
        {
            if (distance_value <= distance_threshold)
                it->second = distance_value;
            else
                this->erase(it);
        }
    }

    value_type get_distance_value(size_t index)
    {
        assert(nseq != 0);

        auto it = this->find(index);

        if (it != this->end())
            return it->second;
        else
            return dummy;
    }
};

inline auto length(map_distance_matrix const & m)
{
    return m.nseq * m.nseq;
}
