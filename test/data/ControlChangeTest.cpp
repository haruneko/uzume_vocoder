// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/ControlChange.hpp"

class ControlChangeTest : public ::testing::Test {
public:
    virtual void setUp() { }
};

using namespace uzume::vocoder;

TEST_F(ControlChangeTest, addedValueShouldBeSortedAndReturnValuesProperly) {
    ControlChange sut(1.0);
    sut.add({0.25, 1.5});
    sut.add({1.0, 0.0});
    sut.add({0.5, 1.0});
    sut.add({0.75, 0.5});

    EXPECT_DOUBLE_EQ(sut.at(0.0), 1.0);
    EXPECT_DOUBLE_EQ(sut.at(0.25), 1.5);
    EXPECT_DOUBLE_EQ(sut.at(0.5), 1.0);
    EXPECT_DOUBLE_EQ(sut.at(0.75), 0.5);
    EXPECT_DOUBLE_EQ(sut.at(1.0), 0.0);
    EXPECT_DOUBLE_EQ(sut.at(0.125), 1.25);
    EXPECT_DOUBLE_EQ(sut.at(0.375), 1.25);
    EXPECT_DOUBLE_EQ(sut.at(0.625), 0.75);
    EXPECT_DOUBLE_EQ(sut.at(0.875), 0.25);
}
