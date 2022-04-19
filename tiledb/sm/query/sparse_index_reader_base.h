/**
 * @file   sparse_index_reader_base.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2021 TileDB, Inc.
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
 * This file defines class SparseIndexReaderBase.
 */

#ifndef TILEDB_SPARSE_INDEX_READER_BASE_H
#define TILEDB_SPARSE_INDEX_READER_BASE_H

#include <queue>
#include "reader_base.h"
#include "tiledb/common/common.h"
#include "tiledb/common/status.h"
#include "tiledb/sm/array_schema/dimension.h"
#include "tiledb/sm/query/query_condition.h"
#include "tiledb/sm/query/result_cell_slab.h"

namespace tiledb {
namespace sm {

class Array;
class ArraySchema;
class MemoryTracker;
class StorageManager;
class Subarray;

/** Processes read queries. */
class SparseIndexReaderBase : public ReaderBase {
 public:
  /* ********************************* */
  /*          TYPE DEFINITIONS         */
  /* ********************************* */

  /** The state for a read sparse global order query. */
  struct ReadState {
    /** The tile index inside of each fragments. */
    std::vector<std::pair<uint64_t, uint64_t>> frag_tile_idx_;

    /** Is the reader done with the query. */
    bool done_adding_result_tiles_;
  };

  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  /** Constructor. */
  SparseIndexReaderBase(
      stats::Stats* stats,
      shared_ptr<Logger> logger,
      StorageManager* storage_manager,
      Array* array,
      Config& config,
      std::unordered_map<std::string, QueryBuffer>& buffers,
      Subarray& subarray,
      Layout layout,
      QueryCondition& condition);

  /** Destructor. */
  ~SparseIndexReaderBase() = default;

  /* ********************************* */
  /*          PUBLIC METHODS           */
  /* ********************************* */

  /**
   * Returns the current read state.
   *
   * @return pointer to the read state.
   */
  const ReadState* read_state() const;

  /**
   * Returns the current read state.
   *
   * @return pointer to the read state.
   */
  ReadState* read_state();

  /**
   * Initializes the reader.
   *
   * @return Status.
   */
  Status init();

  /**
   * Resize the output buffers to the correct size after copying.
   *
   * @param cells_copied Number of cells copied.
   *
   * @return Status.
   */
  Status resize_output_buffers(uint64_t cells_copied);

 protected:
  /* ********************************* */
  /*       PROTECTED ATTRIBUTES        */
  /* ********************************* */

  /** Read state. */
  ReadState read_state_;

  /** Have we loaded all thiles for this fragment. */
  std::vector<uint8_t> all_tiles_loaded_;

  /** Dimension names. */
  std::vector<std::string> dim_names_;

  /** Are dimensions var sized. */
  std::vector<bool> is_dim_var_size_;

  /** Have ve loaded the initial data. */
  bool initial_data_loaded_;

  /** Total memory budget. */
  uint64_t memory_budget_;

  /** Mutex protecting memory budget variables. */
  std::mutex mem_budget_mtx_;

  /** Memory tracker object for the array. */
  MemoryTracker* array_memory_tracker_;

  /** Memory used for coordinates tiles. */
  uint64_t memory_used_for_coords_total_;

  /** How much of the memory budget is reserved for coords. */
  double memory_budget_ratio_coords_;

  /** How much of the memory budget is reserved for array data. */
  double memory_budget_ratio_array_data_;

  /** Are we in elements mode. */
  bool elements_mode_;

  /** Names of dim/attr loaded for query condition. */
  std::vector<std::string> qc_loaded_names_;

  /* Are the users buffers full. */
  bool buffers_full_;

  /* ********************************* */
  /*         PROTECTED METHODS         */
  /* ********************************* */

  /**
   * Return how many cells were copied to the users buffers so far.
   *
   * @param names Attribute/dimensions to compute for.
   *
   * @return Number of cells copied.
   */
  uint64_t cells_copied(const std::vector<std::string>& names);

  /**
   * Load tile offsets.
   *
   * @return Status.
   */
  Status load_offsets();

  /**
   * Read and unfilter coord tiles.
   *
   * @param include_coords Include coordinates or not.
   * @param result_tiles The result tiles to process.
   *
   * @return Status.
   */
  Status read_and_unfilter_coords(
      bool include_coords, const std::vector<ResultTile*>& result_tiles);

  /**
   * Read and unfilter as many attributes as can fit in the memory budget and
   * return the names loaded in 'names_to_copy'. Also keep the 'buffer_idx'
   * updated to keep track of progress.
   *
   * @param memory_budget Memory budget allowed for this operation.
   * @param names Attribute/dimensions to compute for.
   * @param mem_usage_per_attr Computed per attribute memory usage.
   * @param buffer_idx Stores/return the current buffer index in process.
   * @param result_tiles Result tiles to process.
   *
   * @return Status, index_to_copy.
   */
  tuple<Status, optional<std::vector<uint64_t>>> read_and_unfilter_attributes(
      const uint64_t memory_budget,
      const std::vector<std::string>& names,
      const std::vector<uint64_t>& mem_usage_per_attr,
      uint64_t* buffer_idx,
      std::vector<ResultTile*>& result_tiles);

  /**
   * Adds an extra offset in the end of the offsets buffer indicating the
   * returned data size if an attribute is var-sized.
   *
   * @return Status.
   */
  Status add_extra_offset();
};

}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_SPARSE_INDEX_READER_BASE_H
