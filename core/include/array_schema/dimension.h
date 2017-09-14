/**
 * @file   dimension.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017 TileDB, Inc.
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
 * This file defines class Dimension.
 */

#ifndef TILEDB_DIMENSION_H
#define TILEDB_DIMENSION_H

#include <string>

#include "buffer.h"
#include "compressor.h"
#include "datatype.h"
#include "status.h"

namespace tiledb {

/** Manipulates a TileDB dimension. */
class Dimension {
 public:
  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  /** Constructor. */
  Dimension();

    /**
   * Constructor.
   *
   * @param name The name of the dimension.
   * @param type The type of the dimension.
   * @param domain The domain of the dimension.
   * @param tile_extent The tile extent of the dimension.
   */
  Dimension(
      const char* name,
      Datatype type,
      const void* domain,
      const void* tile_extent);

  /**
   * Constructor. It clones the input.
   *
   * @param dim The dimension to clone.
   */
  explicit Dimension(const Dimension* dim);

  /** Destructor. */
  ~Dimension();

  /* ********************************* */
  /*                API                */
  /* ********************************* */

  /** Returns the compressor. */
  Compressor compressor() const;

  /** Returns the compression level. */
  int compression_level() const;

  /**
   * Populates the object members from the data in the input binary buffer.
   *
   * @param buff The buffer to deserialize from.
   * @return Status
   */
  Status deserialize(ConstBuffer* buff);

  /** Returns the domain. */
  void* domain() const;

  /** Dumps the dimension contents in ASCII form in the selected output. */
  void dump(FILE* out) const;

  /** Returns the dimension name. */
  const std::string& name() const;

  /**
   * Serializes the object members into a binary buffer.
   *
   * @param buff The buffer to serialize the data into.
   * @return Status
   */
  Status serialize(Buffer* buff);

  /** Sets the dimension compressor. */
  void set_compressor(Compressor compressor);

  /** Sets the dimension compression level. */
  void set_compression_level(int compression_level);

  /** Returns the tile extent. */
  void* tile_extent() const;

  /** Returns the dimension type. */
  Datatype type() const;

 private:
  /* ********************************* */
  /*         PRIVATE ATTRIBUTES        */
  /* ********************************* */

  /** The dimension compressor. */
  Compressor compressor_;

  /** The dimension compression level. */
  int compression_level_;

  /** The dimension domain. */
  void* domain_;

  /** The dimension name. */
  std::string name_;

  /** The tile extent of the dimension. */
  void* tile_extent_;

  /** The dimension type. */
  Datatype type_;
};

}  // namespace tiledb

#endif  // TILEDB_DIMENSION_H
