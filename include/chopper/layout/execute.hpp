#pragma once

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>

namespace chopper::layout
{

int execute(chopper::configuration & config, std::vector<std::string> const & filenames, chopper::data_store & data);

} // namespace chopper::layout
