/**
 * @file serializable_stats.h
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
 * This file declares the class SerializableStats.
 */

#ifndef TILEDB_SERIALIZABLE_STATS_H
#define TILEDB_SERIALIZABLE_STATS_H

#include "tiledb/sm/serialization/tiledb-rest.h"

#include <string>
#include <unordered_map>

namespace tiledb {
namespace sm {

namespace stats {
class Stats;
}

namespace serialization {

// Old function to get stats from capnp. Should be removed once all objects
// have their proper constructors taking in a created stats pointer rather
// than creating their own.
void stats_from_capnp(
    const capnp::Stats::Reader& stats_reader, stats::Stats* stats);

class SerializableStats {
 public:
  /* ****************************** */
  /*   CONSTRUCTORS & DESTRUCTORS   */
  /* ****************************** */

  /**
   * Constructor, to be called from the Stats object.
   *
   * @param timers Reference to the timers.
   * @param counters Reference to the counters.
   */
  SerializableStats(
      const std::unordered_map<std::string, double>& timers,
      const std::unordered_map<std::string, uint64_t>& counters)
      : timers_(timers)
      , counters_(counters) {
  }

  /* ****************************** */
  /*              API               */
  /* ****************************** */

  /**
   * Serializes the object to Json.
   *
   * @return Json string.
   */
  std::string to_json();

  /**
   * Serializes the object to Capnp.
   *
   * @param builder The Stats Capnp builder.
   */
  void to_capnp(capnp::Stats::Builder& builder);

  /**
   * Deserializes from a Capnp reader and create a Stats object.
   *
   * @param parent_stats Pointer to the parent stats to create a child from.
   * @param reader Capnp reader to deserialize from.
   */
  static stats::Stats* from_capnp(
      stats::Stats* parent_stats, const capnp::Stats::Reader& reader);

 private:
  /* ****************************** */
  /*       PRIVATE ATTRIBUTES       */
  /* ****************************** */

  /** Const reference to the timers. */
  const std::unordered_map<std::string, double>& timers_;

  /** Const reference to the counters. */
  const std::unordered_map<std::string, uint64_t>& counters_;
};

}  // namespace serialization
}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_SERIALIZABLE_STATS_H