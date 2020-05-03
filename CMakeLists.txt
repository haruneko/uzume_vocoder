# Copyright 2020 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.13)

project(uzume_dsp CXX)

add_library(uzume_dsp)
target_compile_features(uzume_dsp
        PRIVATE
        cxx_std_11)
target_compile_options(uzume_dsp
        PRIVATE
        -Wall -Wextra -Wno-long-long -pedantic)
target_sources(uzume_dsp
    PRIVATE
        fft.cpp
        GaussianNoiseGenerator.cpp
        ImpulseResponse.cpp
        Spectrogram.cpp
        Spectrum.cpp
        SynthesizeImpulseResponse.cpp
        SynthesizePhrase.cpp)
target_include_directories(uzume_dsp
    INTERFACE
        ${PROJECT_SOURCE_DIR})
