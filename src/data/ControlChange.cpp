// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include <cstdio>
#include "ControlChange.hpp"

using namespace uzume::vocoder;

ControlChange::ControlChange(const ControlChange &other) : points(other.points) {
}

ControlChange::ControlChange(double initialValue) : points() {
    clear(initialValue);
}

double ControlChange::at(double position) const {
    auto i = points.begin();
    for(; i != points.end(); i++) {
        if(position < i->position) {
            break;
        }
    }
    if(i == points.end()) {
        auto p = i;
        p--;
        return p->ratio;
    } else if(i == points.begin()) {
        return i->ratio;
    } else {
        auto p = i;
        p--;
        double interop = (position - p->position) / (i->position - p->position);
        return p->ratio * (1.0 - interop) + i->ratio * interop;
    }
}

void ControlChange::add(const ControlPoint &p) {
    for(auto i = points.begin(); i != points.end(); i++) {
        if(p.position < i->position) {
            points.insert(i, p);
            return;
        }
    }
    points.push_back(p);
}

void ControlChange::clear(double initialValue) {
    points.clear();
    points.push_back({0.0, initialValue});
}
