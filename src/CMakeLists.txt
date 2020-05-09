# Copyright 2020 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.13)

project(uzume_dsp VERSION 0.0.1 LANGUAGES CXX)

add_library(uzume_dsp STATIC)
set_property(TARGET uzume_dsp PROPERTY CXX_STANDARD 11)
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
set_property(TARGET
        PROPERTY PUBLIC_HEADER
        GaussianNoiseGenerator.hpp
        ImpulseResponse.hpp
        Spectrogram.hpp
        Spectrum.hpp
        SynthesizeImpulseResponse.hpp
        SynthesizePhrase.hpp)
target_include_directories(uzume_dsp
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)
install(TARGETS uzume_dsp
        EXPORT libuzume_dsp
        ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
        PUBLIC_HEADER DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)
