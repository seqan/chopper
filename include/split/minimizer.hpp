#pragma once

#include <seqan/index.h>

struct minimizer
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Attention: all operations on a solely default constructed decorator,
    //!       except assigning a new range, are UB.
    constexpr minimizer() = default;
    //!\brief Copy constructor.
    constexpr minimizer(minimizer const &) = default;
    //!\brief Copy construction via assignment.
    constexpr minimizer & operator=(minimizer const &) = default;
    //!\brief Move constructor.
    constexpr minimizer(minimizer && rhs) = default;
    //!\brief Move assignment.
    constexpr minimizer & operator=(minimizer && rhs) = default;
    //!\brief Use default deconstructor.
    ~minimizer() = default;

    //!\brief Copy constructor from uint64_t.
    constexpr minimizer(uint64_t const v) : value{v} {};
    //!\brief Copy construction via assignment from uint64_t.
    constexpr minimizer & operator=(uint64_t const v) { value = v; return *this; };
    //!\brief Move constructor from uint64_t.
    constexpr minimizer(uint64_t && v) : value{v} {};
    //!\brief Move assignment from uint64_t.
    constexpr minimizer & operator=(uint64_t && v) { value = v; return *this; };

    constexpr minimizer(uint64_t const v, uint64_t const p) : value{v}, position{p} {};

    operator uint64_t() const
    {
        return value;
    }

    uint64_t value{};
    uint64_t position{};

    constexpr friend bool operator==(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value == rhs.value;
    }

    constexpr friend bool operator!=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value != rhs.value;
    }

    constexpr friend bool operator<(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value < rhs.value;
    }
    constexpr friend bool operator<=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value <= rhs.value;
    }
    constexpr friend bool operator>(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value > rhs.value;
    }
    constexpr friend bool operator>=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value >= rhs.value;
    }
};

struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    // Shape for forward hashes
    seqan::Shape<seqan::Dna, seqan::SimpleShape> kmerShape;
    // Shape for hashes on reverse complement
    seqan::Shape<seqan::Dna, seqan::SimpleShape> revCompShape;
    // k-mer size
    uint16_t k{20};
    // window size
    uint32_t w{20};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    template<typename TIt>
    inline void revHashInit(TIt it)
    {
        seqan::hashInit(revCompShape, it);
    }

    template<typename TIt>
    inline auto revHashNext(TIt it)
    {
        return seqan::hashNext(revCompShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline void resize(uint16_t newKmerSize, uint32_t neww)
    {
        k = newKmerSize;
        w = neww;
        seqan::resize(kmerShape, k);
        seqan::resize(revCompShape, k);
    }

    seqan::String<minimizer> getMinimizer(seqan::String<seqan::Dna> const & text)
    {
        if (k > seqan::length(text))
            return seqan::String<minimizer>{};

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;

        seqan::String<minimizer> kmerHashes{};
        reserve(kmerHashes, possible); // maybe rather reserve to expected?

        // Stores hash, begin and end for all k-mers in the window
        std::deque<uint64_t> windowValues;

        auto it = begin(text);
        hashInit(it);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;

            windowValues.push_back(kmerHash);

            ++it;
        }

        auto less_or_equal_compare = [] (auto const & a, auto const & b) { return a <= b; };
        auto min = std::min_element(std::begin(windowValues), std::end(windowValues), less_or_equal_compare);
        seqan::appendValue(kmerHashes, minimizer{*min, static_cast<uint64_t>(std::distance(std::begin(windowValues), min))});
        // appendValue(kmerHashPoss, );

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        uint64_t current_pos{kmerHashes[0].position};
        for (uint64_t pos = 1; pos < possible; ++pos)
        {
            windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            windowValues.push_back(kmerHash);
            ++it;

            if (kmerHash < back(kmerHashes).value)
            {
                current_pos = pos + windowValues.size() - 1; // at end
                seqan::appendValue(kmerHashes, minimizer{kmerHash, current_pos});

            }
            else if (current_pos == pos - 1) // minimum was at the beginning and went out of scope
            {
                min = std::min_element(std::begin(windowValues), std::end(windowValues), less_or_equal_compare);

                if (current_pos != pos + std::distance(std::begin(windowValues), min))
                {
                    current_pos = pos + std::distance(std::begin(windowValues), min);
                    seqan::appendValue(kmerHashes, minimizer{*min, current_pos});
                }
            }
        }

        return kmerHashes;
    }
};
