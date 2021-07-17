// Copyright 2021 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/LinearTimeAxisMap.hpp"

class LinearTimeAxisMapTest : public ::testing::Test {
public:
    virtual void setUp() { }
};

using namespace uzume::vocoder;

TEST_F(LinearTimeAxisMapTest, atShouldRetrunItsRelativeTimePosition) {
    LinearTimeAxisMap sut(100.0, 200.0, 200.0);

    EXPECT_DOUBLE_EQ(sut.at(-50.0), 75.0);
    EXPECT_DOUBLE_EQ(sut.at(0.0), 100.0);
    EXPECT_DOUBLE_EQ(sut.at(50.0), 125.0);
    EXPECT_DOUBLE_EQ(sut.at(100.0), 150.0);
    EXPECT_DOUBLE_EQ(sut.at(150.0), 175.0);
    EXPECT_DOUBLE_EQ(sut.at(200.0), 200.0);
    EXPECT_DOUBLE_EQ(sut.at(250.0), 225.0);
}
