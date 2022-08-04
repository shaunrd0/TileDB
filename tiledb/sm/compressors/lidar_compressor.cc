/**
 * @file   lidar_compressor.cc
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
 * This file implements the lidar compressor class.
 */

#include "tiledb/sm/compressors/lidar_compressor.h"
#include "tiledb/common/logger.h"
#include "tiledb/sm/buffer/buffer.h"
#include "tiledb/sm/enums/datatype.h"
#include "tiledb/sm/filter/filter_buffer.h"
#include "tiledb/sm/filter/filter_storage.h"
#include "tiledb/sm/compressors/gzip_compressor.h"
#include "tiledb/sm/tile/tile.h"
#include "tiledb/sm/filter/xor_filter.h"

#include <cmath>
#include <vector>
#include <utility>

using namespace tiledb::common;

namespace tiledb {
namespace sm {

template<typename W>
Status Lidar::compress(
    Datatype type, int level, ConstBuffer* input_buffer, Buffer* output_buffer) {
  assert(sizeof(W) == 4 || sizeof(W) == 8);
  assert((input_buffer->size() & sizeof(W)) == 0);
  std::vector<std::pair<W, size_t>> vals;
  for (size_t i = 0; i < input_buffer->size()/sizeof(W); ++i) {
    W val = input_buffer->value<W>(i * sizeof(W));
    vals.push_back(std::make_pair(val, i));
  }
  // Sort values
  /// TODO: #include "tiledb/sm/misc/parallel_functions.h" parallel sort
  std::sort(vals.begin(), vals.end());

  std::vector<W> num_vals;
  for (const auto &elem : vals) {
    num_vals.push_back(elem.first);
  }

  // Apply XOR filter.
  // convert std::vector to filter buffer?? 
  // examples of using xor filter in buffer
  Datatype int_type = type == Datatype::FLOAT32 ? Datatype::INT32 : Datatype::INT64;
  Tile tile;
    tile.init_unfiltered(
        constants::format_version, int_type, 0, 1, 0);
  FilterBuffer input;
  input.init(num_vals.data(), num_vals.size() * sizeof(W));
  FilterBuffer output;
  FilterBuffer input_metadata;
  FilterBuffer output_metadata;
  xor_filter_.run_forward(tile, nullptr, &input, &input_metadata, &output, &output_metadata);

  // Apply GZIP compressor.
  assert(output.num_buffers() == 1);
  GZip::compress(level, &output.buffers()[0], output_buffer);
  return Status::Ok();
}

Status Lidar::compress(
    Datatype type, int level, ConstBuffer* input_buffer, Buffer* output_buffer) {
  switch (type) {
    case Datatype::FLOAT32: {
      return compress<int32_t>(type, level, input_buffer, output_buffer);
    }
    case Datatype::FLOAT64: {
      return compress<int64_t>(type, level, input_buffer, output_buffer);
    }
    default: {
      return Status_CompressionError("Lidar::compress: attribute type is not a floating point type.");
    }
  }
}

Status Lidar::compress(Datatype type, ConstBuffer* input_buffer, Buffer* output_buffer) {
    return compress(type, default_level_, input_buffer, output_buffer);
}

Status Lidar::decompress(
    Datatype type, ConstBuffer* input_buffer, PreallocatedBuffer* output_buffer) {
  (void)type;
  (void)input_buffer;
  (void)output_buffer;
  return Status::Ok();
}

uint64_t Lidar::overhead(uint64_t nbytes) {
  (void)nbytes;
  return 0;
}

}  // namespace sm
}  // namespace tiledb
