// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef __DSP_IMPULSE_RESPONSE_HPP__
#define __DSP_IMPULSE_RESPONSE_HPP__

/**
 * ImpulseResponse buffers an impulse response.
 */
class ImpulseResponse final {
public:
    ImpulseResponse(unsigned int length);

    ~ImpulseResponse();

    const double *data() const { return raw; }

    double *data() { return raw; }

    int length() const { return this->_length; }

private:
    ImpulseResponse() {}

    double *raw;
    unsigned int _length;
};

#endif
