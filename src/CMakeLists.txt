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
        world/AnalyzeAperiodicityWithD4C.cpp
        world/AnalyzePeriodicityWithCheapTrick.cpp
        world/audioio.cpp
        world/constant.cpp
        data/Contour.cpp
        world/EstimateF0WithDIO.cpp
        world/fft.cpp
        world/GaussianNoiseGenerator.cpp
        data/ImpulseResponse.cpp
        data/InstantWaveform.cpp
        spectrogram/NaiveSpectrogram.cpp
        data/Spectrum.cpp
        world/SynthesizeImpulseResponseWithWORLD.cpp
        SynthesizeSegment.cpp
        world/SynthesizeSegmentWithWORLD.cpp
        world/util.cpp
        data/Waveform.cpp
        spectrogram/WaveformSpectrogram.cpp)
set_property(TARGET uzume_dsp
        PROPERTY PUBLIC_HEADER
        AnalyzeAperiodicity.hpp
        world/AnalyzeAperiodicityWithD4C.hpp
        AnalyzePeriodicity.hpp
        world/AnalyzePeriodicityWithCheapTrick.hpp
        data/Contour.hpp
        EstimateF0.hpp
        world/EstimateF0WithDIO.hpp
        world/fft.hpp
        world/GaussianNoiseGenerator.hpp
        data/ImpulseResponse.hpp
        data/InstantWaveform.hpp
        spectrogram/NaiveSpectrogram.hpp
        Spectrogram.hpp
        data/Spectrum.hpp
        SynthesizeImpulseResponse.hpp
        world/SynthesizeImpulseResponseWithWORLD.hpp
        SynthesizeSegment.hpp
        world/SynthesizeSegmentWithWORLD.hpp
        world/util.hpp
        data/Waveform.hpp
        spectrogram/WaveformSpectrogram.hpp)
target_include_directories(uzume_dsp
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>)
install(TARGETS uzume_dsp
        EXPORT libuzume_dsp
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/dsp)
