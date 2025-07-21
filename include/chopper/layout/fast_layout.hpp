#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/layout/layout.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>


namespace chopper::layout
{

void fast_layout(chopper::configuration const & config,
                 std::vector<size_t> const & positions,
                 std::vector<size_t> const & cardinalities,
                 std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                 std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                 seqan::hibf::layout::layout & hibf_layout);

} // namespace chopper::layout
