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
target_include_directories(uzume_dsp
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)
install(TARGETS uzume_dsp
        EXPORT uzume_dsp-config)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/uzume_dsp-config-version.cmake
        COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/uzume_dsp-config-version.cmake
        DESTINATION lib/cmake/uzume_dsp)

install(EXPORT uzume_dsp-config
        NAMESPACE uzume::
        DESTINATION lib/cmake/uzume_dsp)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
        DESTINATION include)
