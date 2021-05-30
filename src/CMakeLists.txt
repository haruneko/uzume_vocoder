# Copyright 2020 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.13)

project(uzume_vocoder VERSION 0.0.1 LANGUAGES CXX)

add_library(uzume_vocoder STATIC)
set_property(TARGET uzume_vocoder PROPERTY CXX_STANDARD 11)
target_compile_features(uzume_vocoder
        PRIVATE
        cxx_std_11)
target_compile_options(uzume_vocoder
        PRIVATE
        -Wall -Wextra -Wno-long-long -pedantic)
target_sources(uzume_vocoder
        PRIVATE
        data/Contour.cpp
        data/ImpulseResponse.cpp
        data/InstantWaveform.cpp
        data/Spectrum.cpp
        data/Waveform.cpp
        spectrogram/NaiveSpectrogram.cpp
        spectrogram/StretchedPartialSpectrogram.cpp
        spectrogram/WaveformSpectrogram.cpp
        world/AnalyzeAperiodicityWithD4C.cpp
        world/AnalyzePeriodicityWithCheapTrick.cpp
        world/audioio.cpp
        world/constant.cpp
        world/EstimateF0WithDIO.cpp
        world/fft.cpp
        world/GaussianNoiseGenerator.cpp
        world/SynthesizeImpulseResponseWithWORLD.cpp
        world/SynthesizeSegmentWithWORLD.cpp
        world/util.cpp
        SynthesizeSegment.cpp)
set_property(TARGET uzume_vocoder
        PROPERTY PUBLIC_HEADER
        AnalyzeAperiodicity.hpp
        AnalyzePeriodicity.hpp
        EstimateF0.hpp
        Spectrogram.hpp
        SynthesizeImpulseResponse.hpp
        SynthesizeSegment.hpp
        data/Contour.hpp
        data/ImpulseResponse.hpp
        data/InstantWaveform.hpp
        data/Spectrum.hpp
        data/Waveform.hpp
        spectrogram/NaiveSpectrogram.hpp
        spectrogram/StretchedPartialSpectrogram.hpp
        spectrogram/WaveformSpectrogram.hpp
        world/AnalyzeAperiodicityWithD4C.hpp
        world/AnalyzePeriodicityWithCheapTrick.hpp
        world/EstimateF0WithDIO.hpp
        world/fft.hpp
        world/GaussianNoiseGenerator.hpp
        world/SynthesizeImpulseResponseWithWORLD.hpp
        world/SynthesizeSegmentWithWORLD.hpp
        world/util.hpp)
target_include_directories(uzume_vocoder
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)
install(TARGETS uzume_vocoder
        EXPORT libuzume_vocoder
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(DIRECTORY . DESTINATION include/uzume_vocoder)
