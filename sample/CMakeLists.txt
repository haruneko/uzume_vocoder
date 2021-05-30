# Copyright 2021 Hal@shurabaP.  All rights reserved.
# Use of this source code is governed by a MIT style
# license that can be found in the LICENSE file.
cmake_minimum_required(VERSION 3.13)

project(uzume_vocoder_sample VERSION 0.0.1 LANGUAGES CXX)

add_executable(uzume_vocoder_sample)

include(ExternalProject)

ExternalProject_add(uzume_vocoder
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/uzume_vocoder
        GIT_REPOSITORY https://github.com/haruneko/uzume_vocoder
        GIT_TAG master
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
        )

add_dependencies(uzume_vocoder_sample uzume_vocoder)
link_directories(uzume_vocoder_sample PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_include_directories(uzume_vocoder_sample PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include/uzume_vocoder)
target_link_libraries(uzume_vocoder_sample PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib/libuzume_vocoder.a)

target_sources(uzume_vocoder_sample
        PRIVATE
        main.cpp)
