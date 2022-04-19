/**
 * @file   sparse_index_reader_base.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2022 TileDB, Inc.
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
 * This file implements class SparseIndexReaderBase.
 */

#include "tiledb/sm/query/sparse_index_reader_base.h"
#include "tiledb/common/logger.h"
#include "tiledb/common/memory_tracker.h"
#include "tiledb/sm/array/array.h"
#include "tiledb/sm/array_schema/array_schema.h"
#include "tiledb/sm/filesystem/vfs.h"
#include "tiledb/sm/fragment/fragment_metadata.h"
#include "tiledb/sm/misc/parallel_functions.h"
#include "tiledb/sm/misc/resource_pool.h"
#include "tiledb/sm/query/iquery_strategy.h"
#include "tiledb/sm/query/query_buffer.h"
#include "tiledb/sm/query/query_macros.h"
#include "tiledb/sm/query/strategy_base.h"
#include "tiledb/sm/subarray/subarray.h"

namespace tiledb {
namespace sm {

/* ****************************** */
/*          CONSTRUCTORS          */
/* ****************************** */

SparseIndexReaderBase::SparseIndexReaderBase(
    stats::Stats* stats,
    shared_ptr<Logger> logger,
    StorageManager* storage_manager,
    Array* array,
    Config& config,
    std::unordered_map<std::string, QueryBuffer>& buffers,
    Subarray& subarray,
    Layout layout,
    QueryCondition& condition)
    : ReaderBase(
          stats,
          logger,
          storage_manager,
          array,
          config,
          buffers,
          subarray,
          layout,
          condition)
    , initial_data_loaded_(false)
    , memory_budget_(0)
    , array_memory_tracker_(array->memory_tracker())
    , memory_used_for_coords_total_(0)
    , memory_budget_ratio_coords_(0.5)
    , memory_budget_ratio_array_data_(0.1)
    , buffers_full_(false) {
  read_state_.done_adding_result_tiles_ = false;
}

/* ****************************** */
/*        PROTECTED METHODS       */
/* ****************************** */

const typename SparseIndexReaderBase::ReadState*
SparseIndexReaderBase::read_state() const {
  return &read_state_;
}

typename SparseIndexReaderBase::ReadState* SparseIndexReaderBase::read_state() {
  return &read_state_;
}

Status SparseIndexReaderBase::init() {
  // Sanity checks
  if (storage_manager_ == nullptr)
    return logger_->status(Status_ReaderError(
        "Cannot initialize sparse global order reader; Storage manager not "
        "set"));
  if (buffers_.empty())
    return logger_->status(Status_ReaderError(
        "Cannot initialize sparse global order reader; Buffers not set"));

  // Check subarray
  RETURN_NOT_OK(check_subarray());

  // Load offset configuration options.
  bool found = false;
  offsets_format_mode_ = config_.get("sm.var_offsets.mode", &found);
  assert(found);
  if (offsets_format_mode_ != "bytes" && offsets_format_mode_ != "elements") {
    return logger_->status(
        Status_ReaderError("Cannot initialize reader; Unsupported offsets "
                           "format in configuration"));
  }
  elements_mode_ = offsets_format_mode_ == "elements";

  RETURN_NOT_OK(config_.get<bool>(
      "sm.var_offsets.extra_element", &offsets_extra_element_, &found));
  assert(found);
  RETURN_NOT_OK(config_.get<uint32_t>(
      "sm.var_offsets.bitsize", &offsets_bitsize_, &found));
  if (offsets_bitsize_ != 32 && offsets_bitsize_ != 64) {
    return logger_->status(
        Status_ReaderError("Cannot initialize reader; "
                           "Unsupported offsets bitsize in configuration"));
  }

  // Check the validity buffer sizes.
  RETURN_NOT_OK(check_validity_buffer_sizes());

  return Status::Ok();
}

uint64_t SparseIndexReaderBase::cells_copied(
    const std::vector<std::string>& names) {
  auto& last_name = names.back();
  auto buffer_size = *buffers_[last_name].buffer_size_;
  if (array_schema_.var_size(last_name)) {
    if (buffer_size == 0)
      return 0;
    else
      return buffer_size / (offsets_bitsize_ / 8) - offsets_extra_element_;
  } else {
    return buffer_size / array_schema_.cell_size(last_name);
  }
}

Status SparseIndexReaderBase::load_offsets() {
  // Set a limit to the array memory.
  array_memory_tracker_->set_budget(
      memory_budget_ * memory_budget_ratio_array_data_);

  // Preload zipped coordinate tile offsets. Note that this will
  // ignore fragments with a version >= 5.
  std::vector<std::string> zipped_coords_names = {constants::coords};
  RETURN_CANCEL_OR_ERROR(load_tile_offsets(subarray_, zipped_coords_names));

  // Preload unzipped coordinate tile offsets. Note that this will
  // ignore fragments with a version < 5.
  const auto dim_num = array_schema_.dim_num();
  dim_names_.reserve(dim_num);
  is_dim_var_size_.reserve(dim_num);
  std::vector<std::string> var_size_to_load;
  for (unsigned d = 0; d < dim_num; ++d) {
    dim_names_.emplace_back(array_schema_.dimension(d)->name());
    is_dim_var_size_.emplace_back(array_schema_.var_size(dim_names_[d]));
    if (is_dim_var_size_[d])
      var_size_to_load.emplace_back(dim_names_[d]);
  }
  RETURN_CANCEL_OR_ERROR(load_tile_offsets(subarray_, dim_names_));

  // Compute tile offsets to load and var size to load for attributes.
  std::vector<std::string> attr_tile_offsets_to_load;
  for (auto& it : buffers_) {
    const auto& name = it.first;
    if (array_schema_.is_dim(name))
      continue;

    attr_tile_offsets_to_load.emplace_back(name);

    if (array_schema_.var_size(name))
      var_size_to_load.emplace_back(name);
  }

  // Load tile offsets and var sizes for attributes.
  RETURN_CANCEL_OR_ERROR(load_tile_var_sizes(subarray_, var_size_to_load));
  RETURN_CANCEL_OR_ERROR(
      load_tile_offsets(subarray_, attr_tile_offsets_to_load));

  return Status::Ok();
}

Status SparseIndexReaderBase::read_and_unfilter_coords(
    bool include_coords, const std::vector<ResultTile*>& result_tiles) {
  auto timer_se = stats_->start_timer("read_and_unfilter_coords");

  // Not including coords or no query condition, exit.
  if (!include_coords && condition_.empty())
    return Status::Ok();

  if (subarray_.is_set() || include_coords) {
    // Read and unfilter zipped coordinate tiles. Note that
    // this will ignore fragments with a version >= 5.
    std::vector<std::string> zipped_coords_names = {constants::coords};
    RETURN_CANCEL_OR_ERROR(
        read_coordinate_tiles(zipped_coords_names, result_tiles, true));
    RETURN_CANCEL_OR_ERROR(
        unfilter_tiles(constants::coords, result_tiles, true));

    // Read and unfilter unzipped coordinate tiles. Note that
    // this will ignore fragments with a version < 5.
    RETURN_CANCEL_OR_ERROR(
        read_coordinate_tiles(dim_names_, result_tiles, true));
    for (const auto& dim_name : dim_names_) {
      RETURN_CANCEL_OR_ERROR(unfilter_tiles(dim_name, result_tiles, true));
    }
  }

  if (!condition_.empty()) {
    // Read and unfilter tiles for querty condition.
    RETURN_CANCEL_OR_ERROR(
        read_attribute_tiles(qc_loaded_names_, result_tiles, true));

    for (const auto& name : qc_loaded_names_) {
      RETURN_CANCEL_OR_ERROR(unfilter_tiles(name, result_tiles, true));
    }
  }

  logger_->debug("Done reading and unfiltering coords tiles");
  return Status::Ok();
}

tuple<Status, optional<std::vector<uint64_t>>>
SparseIndexReaderBase::read_and_unfilter_attributes(
    const uint64_t memory_budget,
    const std::vector<std::string>& names,
    const std::vector<uint64_t>& mem_usage_per_attr,
    uint64_t* buffer_idx,
    std::vector<ResultTile*>& result_tiles) {
  auto timer_se = stats_->start_timer("read_and_unfilter_attributes");

  std::vector<std::string> names_to_read;
  std::vector<uint64_t> index_to_copy;
  uint64_t memory_used = 0;
  while (*buffer_idx < names.size()) {
    auto& name = names[*buffer_idx];
    auto attr_mem_usage = mem_usage_per_attr[*buffer_idx];
    if (memory_used + attr_mem_usage < memory_budget) {
      memory_used += attr_mem_usage;

      // We only read attributes, so dimensions have 0 cost.
      if (attr_mem_usage != 0)
        names_to_read.emplace_back(name);

      index_to_copy.emplace_back(*buffer_idx);
      (*buffer_idx)++;
    } else {
      break;
    }
  }

  // Read and unfilter tiles.
  RETURN_NOT_OK_TUPLE(
      read_attribute_tiles(names_to_read, result_tiles, true), nullopt);

  for (auto& name : names_to_read)
    RETURN_NOT_OK_TUPLE(unfilter_tiles(name, result_tiles, true), nullopt);

  return {Status::Ok(), std::move(index_to_copy)};
}

Status SparseIndexReaderBase::resize_output_buffers(uint64_t cells_copied) {
  // Resize buffers if the result cell slabs was truncated.
  for (auto& it : buffers_) {
    const auto& name = it.first;
    const auto size = *it.second.buffer_size_;
    uint64_t num_cells = 0;

    if (array_schema_.var_size(name)) {
      // Get the current number of cells from the offsets buffer.
      num_cells = size / constants::cell_var_offset_size;

      // Remove an element if the extra element flag is set.
      if (offsets_extra_element_ && num_cells > 0)
        num_cells--;

      // Buffer should be resized.
      if (num_cells > cells_copied) {
        // Offsets buffer is trivial.
        *(it.second.buffer_size_) =
            cells_copied * constants::cell_var_offset_size +
            offsets_extra_element_ * offsets_bytesize();

        // Since the buffer is shrunk, there is an offset for the next element
        // loaded, use it.
        if (offsets_bitsize_ == 64) {
          uint64_t offset_div =
              elements_mode_ ? datatype_size(array_schema_.type(name)) : 1;
          *it.second.buffer_var_size_ =
              ((uint64_t*)it.second.buffer_)[cells_copied] * offset_div;
        } else {
          uint32_t offset_div =
              elements_mode_ ? datatype_size(array_schema_.type(name)) : 1;
          *it.second.buffer_var_size_ =
              ((uint32_t*)it.second.buffer_)[cells_copied] * offset_div;
        }
      }
    } else {
      // Always adjust the size for fixed size attributes.
      auto cell_size = array_schema_.cell_size(name);
      *(it.second.buffer_size_) = cells_copied * cell_size;
    }

    // Always adjust validity vector size, if present.
    if (num_cells > cells_copied) {
      if (it.second.validity_vector_.buffer_size() != nullptr)
        *(it.second.validity_vector_.buffer_size()) =
            num_cells * constants::cell_validity_size;
    }
  }

  return Status::Ok();
}

Status SparseIndexReaderBase::add_extra_offset() {
  for (const auto& it : buffers_) {
    const auto& name = it.first;
    if (!array_schema_.var_size(name))
      continue;

    // Do not apply offset for empty results because we will
    // write backwards and corrupt memory we don't own.
    if (*it.second.buffer_size_ == 0)
      continue;

    auto buffer = static_cast<unsigned char*>(it.second.buffer_);
    if (offsets_format_mode_ == "bytes") {
      memcpy(
          buffer + *it.second.buffer_size_ - offsets_bytesize(),
          it.second.buffer_var_size_,
          offsets_bytesize());
    } else if (offsets_format_mode_ == "elements") {
      auto elements =
          *it.second.buffer_var_size_ / datatype_size(array_schema_.type(name));
      memcpy(
          buffer + *it.second.buffer_size_ - offsets_bytesize(),
          &elements,
          offsets_bytesize());
    } else {
      return logger_->status(Status_ReaderError(
          "Cannot add extra offset to buffer; Unsupported offsets format"));
    }
  }

  return Status::Ok();
}

}  // namespace sm
}  // namespace tiledb
