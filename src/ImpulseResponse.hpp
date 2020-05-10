// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_IMPULSE_RESPONSE_HPP__
#define __DSP_IMPULSE_RESPONSE_HPP__

namespace uzume { namespace dsp {

/**
 * ImpulseResponse buffers an impulse response.
 */
class ImpulseResponse final {
public:
    explicit ImpulseResponse(unsigned int length);

    ~ImpulseResponse();

    const double *data() const { return raw; }

    double *data() { return raw; }

    unsigned int length() const { return _length; }

private:
    ImpulseResponse() = default;

    double *raw;
    unsigned int _length;
};

} }

#endif
