# Copyright 2020 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.9)

project(uzume_dsp CXX)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -s")

file(GLOB DSP_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} *.cpp)

add_library(uzume_dsp ${DSP_SOURCES})
