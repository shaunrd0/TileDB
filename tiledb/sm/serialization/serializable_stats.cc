/**
 * @file serializable_stats.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines the class SerializableStats.
 */

// clang-format off
#ifdef TILEDB_SERIALIZATION
#include <capnp/compat/json.h>
#include <capnp/message.h>
#include <capnp/serialize.h>
#include "tiledb/sm/serialization/capnp_utils.h"
#include "serializable_stats.h"
#endif
// clang-format on

#include "tiledb/common/common.h"

namespace tiledb {
namespace sm {
namespace serialization {

#ifdef TILEDB_SERIALIZATION

void stats_from_capnp(
    const capnp::Stats::Reader& stats_reader, stats::Stats* stats) {
  if (stats_reader.hasCounters()) {
    auto counters = stats->counters();
    auto counters_reader = stats_reader.getCounters();
    for (const auto entry : counters_reader.getEntries()) {
      (*counters)[std::string(entry.getKey().cStr())] = entry.getValue();
    }
  }

  if (stats_reader.hasTimers()) {
    auto timers = stats->timers();
    auto timers_reader = stats_reader.getTimers();
    for (const auto entry : timers_reader.getEntries()) {
      (*timers)[std::string(entry.getKey().cStr())] = entry.getValue();
    }
  }
}

std::string SerializableStats::to_json() {
  ::capnp::MallocMessageBuilder message;
  auto builder = message.initRoot<capnp::Stats>();

  to_capnp(builder);
  ::capnp::JsonCodec json;
  kj::String capnp_json = json.encode(builder);

  return std::string(capnp_json.cStr(), capnp_json.size());
}

void SerializableStats::to_capnp(capnp::Stats::Builder& stats_builder) {
  // Build counters
  if (!counters_.empty()) {
    auto counters_builder = stats_builder.initCounters();
    auto entries_builder = counters_builder.initEntries(counters_.size());
    uint64_t index = 0;
    for (const auto& entry : counters_) {
      entries_builder[index].setKey(entry.first);
      entries_builder[index].setValue(entry.second);
      ++index;
    }
  }

  // Build timers
  if (!timers_.empty()) {
    auto timers_builder = stats_builder.initTimers();
    auto entries_builder = timers_builder.initEntries(timers_.size());
    uint64_t index = 0;
    for (const auto& entry : timers_) {
      entries_builder[index].setKey(entry.first);
      entries_builder[index].setValue(entry.second);
      ++index;
    }
  }
}

stats::Stats* SerializableStats::from_capnp(
    stats::Stats* parent_stats, const capnp::Stats::Reader& stats_reader) {
  std::unordered_map<std::string, uint64_t> counters;
  if (stats_reader.hasCounters()) {
    auto counters_reader = stats_reader.getCounters();
    for (const auto entry : counters_reader.getEntries()) {
      counters[std::string(entry.getKey().cStr())] = entry.getValue();
    }
  }

  std::unordered_map<std::string, double> timers;
  if (stats_reader.hasTimers()) {
    auto timers_reader = stats_reader.getTimers();
    for (const auto entry : timers_reader.getEntries()) {
      timers[std::string(entry.getKey().cStr())] = entry.getValue();
    }
  }

  return parent_stats->create_child("Subarray", move(timers), move(counters));
}

#endif

}  // namespace serialization
}  // namespace sm
}  // namespace tiledb
