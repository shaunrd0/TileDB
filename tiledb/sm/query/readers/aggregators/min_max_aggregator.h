/**
 * @file   min_max_aggregator.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 * This file defines class MinMaxAggregator.
 */

#ifndef TILEDB_MIN_MAX_AGGREGATOR_H
#define TILEDB_MIN_MAX_AGGREGATOR_H

#include "tiledb/sm/query/readers/aggregators/aggregate_with_count.h"
#include "tiledb/sm/query/readers/aggregators/iaggregator.h"
#include "tiledb/sm/query/readers/aggregators/min_max.h"
#include "tiledb/sm/query/readers/aggregators/validity_policies.h"

#include <functional>

namespace tiledb::sm {

class QueryBuffer;

/**
 * Comparator aggregator base class to handle partial specialization of some
 * types.
 */
template <typename T>
class ComparatorAggregatorBase {
 protected:
  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  ComparatorAggregatorBase() = delete;
  ComparatorAggregatorBase(const FieldInfo&& field_info) = delete;

  /**
   * Constructor.
   *
   * @param field_info Field info.
   */
  ComparatorAggregatorBase(const FieldInfo& field_info)
      : field_info_(field_info)
      , value_(nullopt)
      , validity_value_(
            field_info_.is_nullable_ ? std::make_optional(0) : nullopt) {
  }

  DISABLE_COPY_AND_COPY_ASSIGN(ComparatorAggregatorBase);
  DISABLE_MOVE_AND_MOVE_ASSIGN(ComparatorAggregatorBase);

  /* ********************************* */
  /*        PROTECTED METHODS          */
  /* ********************************* */

  /**
   * Copy final data to the user buffer.
   *
   * @param output_field_name Name for the output buffer.
   * @param buffers Query buffers.
   */
  void copy_to_user_buffer(
      std::string output_field_name,
      std::unordered_map<std::string, QueryBuffer>& buffers) const;

 protected:
  /* ********************************* */
  /*       PROTECTED ATTRIBUTES        */
  /* ********************************* */

  /** Field information. */
  const FieldInfo field_info_;

  /** Computed min/max. */
  optional<T> value_;

  /** Computed validity value. */
  optional<uint8_t> validity_value_;
};

template <typename T, typename Op>
class ComparatorAggregator : public ComparatorAggregatorBase<T>,
                             public OutputBufferValidator,
                             public IAggregator {
 protected:
  using VALUE_T = typename type_data<T>::value_type;

  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  ComparatorAggregator() = delete;
  ComparatorAggregator(const FieldInfo&& field_info) = delete;

  /**
   * Constructor.
   *
   * @param field_info Field info.
   */
  ComparatorAggregator(const FieldInfo& field_info);

  DISABLE_COPY_AND_COPY_ASSIGN(ComparatorAggregator);
  DISABLE_MOVE_AND_MOVE_ASSIGN(ComparatorAggregator);

  /* ********************************* */
  /*                API                */
  /* ********************************* */

 public:
  /** Returns the field name for the aggregator. */
  std::string field_name() override {
    return ComparatorAggregatorBase<T>::field_info_.name_;
  }

  /** Returns if the aggregation is var sized or not. */
  bool var_sized() override {
    return ComparatorAggregatorBase<T>::field_info_.var_sized_;
  };

  /** Returns if the aggregate needs to be recomputed on overflow. */
  bool need_recompute_on_overflow() override {
    return false;
  }

  /**
   * Validate the result buffer.
   *
   * @param output_field_name Name for the output buffer.
   * @param buffers Query buffers.
   */
  void validate_output_buffer(
      std::string output_field_name,
      std::unordered_map<std::string, QueryBuffer>& buffers) override;

  /**
   * Aggregate data using the aggregator.
   *
   * @param input_data Input data for aggregation.
   */
  void aggregate_data(AggregateBuffer& input_data) override;

  /**
   * Copy final data to the user buffer.
   *
   * @param output_field_name Name for the output buffer.
   * @param buffers Query buffers.
   */
  void copy_to_user_buffer(
      std::string output_field_name,
      std::unordered_map<std::string, QueryBuffer>& buffers) override;

  /** Returns the TileDB datatype of the output field for the aggregate. */
  Datatype output_datatype() override {
    return ComparatorAggregatorBase<T>::field_info_.type_;
  }

 private:
  /* ********************************* */
  /*         PRIVATE ATTRIBUTES        */
  /* ********************************* */

  /** AggregateWithCount to do summation of AggregateBuffer data. */
  AggregateWithCount<T, VALUE_T, MinMax<Op>, NonNull> aggregate_with_count_;

  /** Mutex protecting `value_`. */
  std::mutex value_mtx_;

  /** Operation function to determine value. */
  Op op_;
};

template <typename T>
class MinAggregator : public ComparatorAggregator<
                          T,
                          std::less<typename type_data<T>::value_type>> {
 public:
  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  MinAggregator() = delete;

  /**
   * Constructor.
   *
   * @param field info Field info.
   */
  MinAggregator(const FieldInfo field_info)
      : ComparatorAggregator<T, std::less<typename type_data<T>::value_type>>(
            field_info){};

  DISABLE_COPY_AND_COPY_ASSIGN(MinAggregator);
  DISABLE_MOVE_AND_MOVE_ASSIGN(MinAggregator);
};

template <typename T>
class MaxAggregator : public ComparatorAggregator<
                          T,
                          std::greater<typename type_data<T>::value_type>> {
 public:
  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  MaxAggregator() = delete;

  /**
   * Constructor.
   *
   * @param field info Field info.
   */
  MaxAggregator(const FieldInfo field_info)
      : ComparatorAggregator<
            T,
            std::greater<typename type_data<T>::value_type>>(field_info){};

  DISABLE_COPY_AND_COPY_ASSIGN(MaxAggregator);
  DISABLE_MOVE_AND_MOVE_ASSIGN(MaxAggregator);
};

}  // namespace tiledb::sm

#endif  // TILEDB_MIN_MAX_AGGREGATOR_H
