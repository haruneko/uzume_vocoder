# Copyright 2020 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.16)

project(uzume_vocoder_test)

set(CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.10.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(uzume_vocoder_test)

target_sources(uzume_vocoder_test
        PRIVATE
        data/LinearTimeAxisMapTest.cpp
        data/WaveformTest.cpp
        spectrogram/ArraySpectrogramAggregatorTest.cpp
        spectrogram/StretchedPartialSpectrogramTest.cpp
        spectrogram/WaveformSpectrogramTest.cpp)

target_link_libraries(uzume_vocoder_test
        PRIVATE
        uzume_vocoder
        gmock
        gtest
        gtest_main)

add_test(NAME uzume_vocoder_test COMMAND uzume_vocoder)
