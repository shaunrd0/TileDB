#
# tiledb/sm/serialization/CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2022 TileDB, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

include(common NO_POLICY_SCOPE)

if (TILEDB_TESTS AND TILEDB_SERIALIZATION)
    add_executable(unit_serialization EXCLUDE_FROM_ALL)

    target_compile_definitions(unit_serialization PRIVATE
      -DTILEDB_TEST_INPUTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../../test/inputs"
    )

    # Sources for tests
    target_sources(unit_serialization PUBLIC test/main.cc test/unit_serializable_objects.cc)

    # The dependencies can't yet be factored into separate object libraries
    target_link_libraries(unit_serialization PUBLIC TILEDB_CORE_OBJECTS)
    target_link_libraries(unit_serialization PUBLIC TILEDB_CORE_OBJECTS_ILIB)

    find_package(Catch_EP REQUIRED)
    target_link_libraries(unit_serialization PUBLIC Catch2::Catch2)

    # We want tests to continue as normal even as the API is changing,
    # so don't warn for deprecations, since they'll be escalated to errors.
    if (NOT MSVC)
      target_compile_options(unit_serialization PRIVATE -Wno-deprecated-declarations)
    endif()

    add_test(
        NAME "unit_serialization"
        COMMAND $<TARGET_FILE:unit_serialization>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
