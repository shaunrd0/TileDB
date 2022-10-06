/**
 * @file bitsort_filter.cc
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
 * This file implements the bitsort filter class.
 */

#include "tiledb/common/status.h"
#include "tiledb/common/types/untyped_datum.h"
#include "tiledb/sm/filter/bitsort_filter.h"
#include "tiledb/common/logger.h"
#include "tiledb/sm/buffer/buffer.h"
#include "tiledb/sm/enums/datatype.h"
#include "tiledb/sm/filter/filter_buffer.h"
#include "tiledb/sm/filter/filter_storage.h"
#include "tiledb/sm/misc/constants.h"
#include "tiledb/sm/tile/tile.h"
#include "tiledb/sm/query/query_buffer.h"
#include "tiledb/sm/query/writers/domain_buffer.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace tiledb::common;

namespace tiledb {
namespace sm {

void BitSortFilter::dump(FILE* out) const {
  if (out == nullptr)
    out = stdout;
  fprintf(out, "BitSortFilter");
}

Status BitSortFilter::run_forward(
    const Tile& tile,
    Tile* const tile_offsets,
    FilterBuffer* input_metadata,
    FilterBuffer* input,
    FilterBuffer* output_metadata,
    FilterBuffer* output) const {
  (void)tile;
  (void)tile_offsets;
  (void)input_metadata;
  (void)input;
  (void)output_metadata;
  (void)output;
  return Status_FilterError("BitSortFilter: Do not call default version of run_forward.");
}

Status BitSortFilter::run_reverse(
      const Tile& tile,
      Tile* const tile_offsets,
      FilterBuffer* input_metadata,
      FilterBuffer* input,
      FilterBuffer* output_metadata,
      FilterBuffer* output,
      const Config& config) const {
  (void)tile;
  (void)tile_offsets;
  (void)input_metadata;
  (void)input;
  (void)output_metadata;
  (void)output;
  (void)config;
  return Status_FilterError("BitSortFilter: Do not call default version of of run_reverse.");
}

Status BitSortFilter::run_forward(
    const Tile& tile,
    std::vector<Tile*> &dim_tiles,
    FilterBuffer* input_metadata,
    FilterBuffer* input,
    FilterBuffer* output_metadata,
    FilterBuffer* output) const {

  // Since run_forward interprets the filter's data as integers, we case on
  // the size of the type and pass in the corresponding integer type into
  // a templated function.
  auto tile_type = tile.type();
  switch (datatype_size(tile_type)) {
    case sizeof(uint8_t): {
      return run_forward<uint8_t>(
          dim_tiles, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint16_t): {
      return run_forward<uint16_t>(
         dim_tiles, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint32_t): {
      return run_forward<uint32_t>(
          dim_tiles, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint64_t): {
      return run_forward<uint64_t>(
          dim_tiles, input_metadata, input, output_metadata, output);
    }
    default: {
      return Status_FilterError(
          "BitSortFilter::run_forward: datatype does not have an appropriate "
          "size");
    }
  } 

  return Status_FilterError("BitSortFilter::run_forward: invalid datatype.");
}

template <typename AttrType>
Status BitSortFilter::run_forward(
    std::vector<Tile*> &dim_tiles,
    FilterBuffer* input_metadata,
    FilterBuffer* input,
    FilterBuffer* output_metadata,
    FilterBuffer* output) const {
  // Output size does not change with this filter.
  RETURN_NOT_OK(output->prepend_buffer(input->size()));
  Buffer* output_buf = output->buffer_ptr(0);
  assert(output_buf != nullptr);

  // Write the metadata.
  auto parts = input->buffers();
  auto num_parts = (uint32_t)parts.size();
  uint32_t metadata_size = sizeof(uint32_t) + num_parts * sizeof(uint32_t);
  RETURN_NOT_OK(output_metadata->append_view(input_metadata));
  RETURN_NOT_OK(output_metadata->prepend_buffer(metadata_size));
  RETURN_NOT_OK(output_metadata->write(&num_parts, sizeof(uint32_t)));

  // Keep track of the starting offsets to write in positions for each part.
  std::vector<uint32_t> offsets;
  uint32_t total_size = 0;
  for (const auto& part : parts) {
    auto part_size = (uint32_t)part.size();
    offsets.push_back(total_size);
    total_size += part_size;
  }

  // Sort all parts.
  std::vector<std::pair<AttrType, uint64_t>> sorted_elements(total_size / sizeof(AttrType));
  for (uint64_t i = 0; i < num_parts; ++i) {
    const auto &part = parts[i];
    auto part_size = static_cast<uint32_t>(part.size());
    RETURN_NOT_OK(output_metadata->write(&part_size, sizeof(uint32_t)));
    RETURN_NOT_OK(sort_part<AttrType>(&part, output_buf, offsets[i], sorted_elements));
  }

  // Since rewrite_dim_tile_forward is only moving around dimension data and 
  // reordering it, the type becomes irrelevant and only the size of the type
  // matters, so the code cases on the size of the type and passes in the
  // corresponding integer type into a templated function.
  for (auto* dim_tile : dim_tiles) {
    Datatype tile_type = dim_tile->type();
    switch (datatype_size(tile_type)) {
      case sizeof(uint8_t): {
        rewrite_dim_tile_forward<AttrType, uint8_t>(sorted_elements, dim_tile);
      } break;
      case sizeof(uint16_t): {
        rewrite_dim_tile_forward<AttrType, uint16_t>(sorted_elements, dim_tile);
      } break;
      case sizeof(uint32_t): {
        rewrite_dim_tile_forward<AttrType, uint32_t>(sorted_elements, dim_tile);
      } break;
      case sizeof(uint64_t): {
        rewrite_dim_tile_forward<AttrType, uint64_t>(sorted_elements, dim_tile);
      } break;
      default: {
        return Status_FilterError(
          "BitSortFilter::sort_part: dimension datatype does not have an appropriate "
          "size");
      }
    }
  }

  return Status::Ok();
}

template <typename AttrType>
Status BitSortFilter::sort_part(const ConstBuffer* input_buffer, Buffer* output_buffer, uint32_t start, std::vector<std::pair<AttrType, uint64_t>> &sorted_elements) const {
  // Read in the data.
  uint32_t s = input_buffer->size();
  assert(s % sizeof(AttrType) == 0);
  uint32_t num_elems_in_part = s / sizeof(AttrType);

  if (num_elems_in_part == 0) {
    return Status::Ok();
  }
  const AttrType* part_array = static_cast<const AttrType*>(input_buffer->data());

  // Create the array to sort by keeping track of elements and positions.
  for (uint32_t i = 0; i < num_elems_in_part; ++i) {
    sorted_elements[start + i] = std::make_pair(part_array[i], start + i);
  }

  // Sort the data.
  std::sort(sorted_elements.begin() + start, sorted_elements.begin() + start + num_elems_in_part);

  // Write in the sorted order to output.
  for (uint32_t j = 0; j < num_elems_in_part; ++j) {
    AttrType value = sorted_elements[start + j].first;
    RETURN_NOT_OK(output_buffer->write(&value, sizeof(AttrType)));

    if (j != num_elems_in_part - 1) {
      output_buffer->advance_offset(sizeof(AttrType));
    }
  }

  return Status::Ok();
}

template <typename AttrType, typename DimType>
Status BitSortFilter::rewrite_dim_tile_forward(const std::vector<std::pair<AttrType, uint64_t>> &elements, Tile *dim_tile) const {
  // Obtain the pointer to the data the code modifies.
  uint64_t elements_size = elements.size();
  std::vector<DimType> tile_data_vec(elements_size);
  DimType *tile_data = static_cast<DimType*>(dim_tile->data());

  // Keep track of the data we should write to the tile.
  for (uint64_t i = 0; i < elements_size; ++i) {
    tile_data_vec[i] = tile_data[elements[i].second];
  }

  // Overwrite the tile.
  RETURN_NOT_OK(dim_tile->write(tile_data_vec.data(), 0, sizeof(DimType) * elements_size));
  FilteredBuffer &filtered_buffer = dim_tile->filtered_buffer();
  filtered_buffer.expand(elements_size * sizeof(DimType));
  memcpy(filtered_buffer.data(), tile_data_vec.data(), sizeof(DimType) * elements_size);

  return Status::Ok();
}

Status BitSortFilter::run_reverse(
    const Tile& tile,
    BitSortFilterMetadataType &pair,
    FilterBuffer* input_metadata,
    FilterBuffer* input,
    FilterBuffer* output_metadata,
    FilterBuffer* output,
    const Config& config) const {
  (void)config;

  // Since run_reverse interprets the filter's data as integers, we case on
  // the size of the type and pass in the corresponding integer type into
  // a templated function.
  auto tile_type = tile.type();
  switch (datatype_size(tile_type)) {
    case sizeof(uint8_t): {
      return run_reverse<uint8_t>(
          pair, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint16_t): {
      return run_reverse<uint16_t>(
          pair, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint32_t): {
      return run_reverse<uint32_t>(
          pair, input_metadata, input, output_metadata, output);
    }
    case sizeof(uint64_t): {
      return run_reverse<uint64_t>(
          pair, input_metadata, input, output_metadata, output);
    }
    default: {
      return Status_FilterError(
          "BitSortFilter::run_reverse: datatype does not have an appropriate "
          "size");
    }
  }

  return Status_FilterError("BitSortFilter::run_reverse: invalid datatype.");
}

template <typename AttrType>
Status BitSortFilter::run_reverse(
    BitSortFilterMetadataType &pair,
    FilterBuffer* input_metadata,
    FilterBuffer* input,
    FilterBuffer* output_metadata,
    FilterBuffer* output) const  {
  // Get number of parts.
  uint32_t num_parts;
  RETURN_NOT_OK(input_metadata->read(&num_parts, sizeof(uint32_t)));

  RETURN_NOT_OK(output->prepend_buffer(input->size()));
  Buffer* output_buf = output->buffer_ptr(0);
  assert(output_buf != nullptr);

  // Determine the positions vector and rewrite the dimension tiles using these positions.
  std::vector<uint64_t> positions;
  rewrite_dim_tiles_reverse(pair, positions);

  // Unsort the attribute data with the positions determined by the dimension global sorting.
  for (uint32_t i = 0; i < num_parts; i++) {
    uint32_t part_size;
    RETURN_NOT_OK(input_metadata->read(&part_size, sizeof(uint32_t)));
    ConstBuffer part(nullptr, 0);
    RETURN_NOT_OK(input->get_const_buffer(part_size, &part));

    RETURN_NOT_OK(unsort_part<AttrType>(positions, &part, output_buf));

    if (output_buf->owns_data()) {
      output_buf->advance_size(part_size);
    }
    output_buf->advance_offset(part_size);
    input->advance_offset(part_size);
  }

  // Output metadata is a view on the input metadata, skipping what was used
  // by this filter.
  auto md_offset = input_metadata->offset();
  RETURN_NOT_OK(output_metadata->append_view(
      input_metadata, md_offset, input_metadata->size() - md_offset));

  return Status::Ok();
}

template <typename AttrType>
Status BitSortFilter::unsort_part(
    const std::vector<uint64_t> &positions, ConstBuffer* input_buffer, Buffer* output_buffer) const {
  // Read in data.
  uint32_t s = input_buffer->size();
  assert(s % sizeof(AttrType) == 0);
  uint32_t num_elems_in_part = s / sizeof(AttrType);

  if (num_elems_in_part == 0) {
    return Status::Ok();
  }

  const AttrType* input_array = static_cast<const AttrType*>(input_buffer->data());
  AttrType* output_array = static_cast<AttrType*>(output_buffer->cur_data());

  // Write in data in original order.
  for (uint32_t i = 0; i < num_elems_in_part; ++i) {
    output_array[i] = input_array[positions[i]];
  }

  return Status::Ok();
}

Status BitSortFilter::rewrite_dim_tiles_reverse(BitSortFilterMetadataType &pair, std::vector<uint64_t> &positions) const {
  std::vector<Tile*> &dim_tiles = pair.first.get();
  auto cmp_fn = pair.second;

  for (uint64_t i = 0; i < dim_tiles[0]->cell_num(); ++i) {
    positions.push_back(i);
  }

  // Sort the dimension data and determine the positions of the dimension tile data.
  std::sort(positions.begin(), positions.end(), cmp_fn);

  // Rewrite the individual tiles with the position vector.
  for (auto *dim_tile : dim_tiles) {
  // Since rewrite_dim_tile_reverse is only moving around dimension data and 
  // reordering it, the type becomes irrelevant and only the size of the type
  // matters, so the code cases on the size of the type and passes in the
  // corresponding integer type into a templated function.
    Datatype tile_type = dim_tile->type();
    switch (datatype_size(tile_type)) {
      case sizeof(uint8_t): {
        rewrite_dim_tile_reverse<uint8_t>(dim_tile, positions);
      } break;
      case sizeof(uint16_t): {
        rewrite_dim_tile_reverse<uint16_t>(dim_tile, positions);
      } break;
      case sizeof(uint32_t): {
        rewrite_dim_tile_reverse<uint32_t>(dim_tile, positions);
      } break;
      case sizeof(uint64_t): {
        rewrite_dim_tile_reverse<uint64_t>(dim_tile, positions);
      } break;
      default: {
        return Status_FilterError(
          "BitSortFilter::rewrite_dim_tiles_reverse: dimension datatype does not have an appropriate "
          "size");
      }
    }
  }

  return Status::Ok();
}

template<typename DimType>
Status BitSortFilter::rewrite_dim_tile_reverse(Tile *dim_tile, std::vector<uint64_t> &positions) const {
  // Obtain the pointer to the data the code modifies.
  uint64_t positions_size = positions.size();
  std::vector<DimType> tile_data_vec(positions_size);
  FilteredBuffer &filtered_buffer = dim_tile->filtered_buffer();
  DimType *tile_data = reinterpret_cast<DimType*>(filtered_buffer.data());

  // Keep track of the data we should write to the tile.
  for (uint64_t i = 0; i < positions_size; ++i) {
    tile_data_vec[i] = tile_data[positions[i]];
  }

  // Overwrite the tile.
  memcpy(dim_tile->data(), tile_data_vec.data(), sizeof(DimType) * positions_size);

  return Status::Ok();
}

BitSortFilter* BitSortFilter::clone_impl() const {
  return tdb_new(BitSortFilter);
}

}  // namespace sm
}  // namespace tiledb