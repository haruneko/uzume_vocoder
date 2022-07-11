// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_CONTROL_CHANGE_HPP
#define UZUME_VOCODER_CONTROL_CHANGE_HPP

#include <list>

namespace uzume { namespace vocoder {

struct ControlPoint final {
    double position;
    double ratio;
};

class ControlChange final {
public:
    ControlChange() = delete;
    explicit ControlChange(const ControlChange &other);
    explicit ControlChange(double initialValue);

    double at(double position) const;
    void add(const ControlPoint &p);
    void clear(double initialValue);
private:
    std::list<ControlPoint> points;
};

} }

#endif //UZUME_VOCODER_CONTROL_CHANGE_HPP
