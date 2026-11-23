#pragma once

#include <unordered_set>

#include "Types/DetectorEvent.hpp"
#include "detail/PairHash.hpp"

/// Hit-related data types
using Cluster = E320Io::ChipEvent::Cluster;
using Hit = std::pair<std::uint16_t, std::uint16_t>;
using HitSet = std::unordered_set<Hit, PairHash>;
