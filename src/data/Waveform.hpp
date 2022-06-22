// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_WAVEFORM_HPP
#define UZUME_VOCODER_WAVEFORM_HPP

#include <vector>

namespace uzume { namespace vocoder {

/**
 * Waveform represents a total waveform.
 */
class Waveform {
public:
    Waveform(unsigned int length, unsigned int samplingFrequency);
    ~Waveform();

    /**
     * @return length of this waveform in milli seconds.
     */
    inline double msLength() const {
        return (double)length / samplingFrequency * 1000.0;
    }

    double *data;
    unsigned int length;
    unsigned int samplingFrequency;

    /**
     * `read` creates Waveform from local a *.wav file.
     * @param filepath is a path of the *.wav file to read.
     * @return nullptr when `read` fails reading otherwise return a Waveform instance.
     */
    static Waveform *read(const char *filepath);

    /**
     * `save` writes a *.wav file that contains waveform.
     * @param filepath to write.
     * @param bits represents wav file bits.
     * @return true when `save` succeeds othewise false.
     */
    bool save(const char *filepath, int bits = 16) const;

    /**
     * `clear` clears waveform data with 0 value.
     */
    void clear();

    /**
     * @param msBegin represents milliseconds of start position.
     * @param msEnd represents milliseconds of end position.
     * @return root mean square between `msBegin` and `msEnd`.
     */
    double rootMeanSquareBetween(double msBegin, double msEnd) const;

    /**
     * @param indexBegin represents the index of start position.
     * @param indexEnd represents the index of end position.
     * @return root mean square between `indexBegin` and `indexEnd`.
     */
    double rootMeanSquareBetween(int indexBegin, int indexEnd) const;

    /**
     * @param msBegin represents milliseconds of start position.
     * @param msEnd represents milliseconds of end position.
     * @return the biggest absolute value between `msBegin` and `msEnd`.
     */
    double maxAbsoluteValueBetween(double msBegin, double msEnd) const;

    /**
     * @param indexBegin represents the index of start position.
     * @param indexEnd represents the index of end position.
     * @return the biggest absolute value between `indexBegin` and `indexEnd`.
     */
    double maxAbsoluteValueBetween(int indexBegin, int indexEnd) const;

    /**
     * @param ms represents milliseconds.
     * @return the corresponding index of waveform array.
     */
    inline int indexAt(double ms) const {
            return (int)(ms / msLength() * length);
    }

    /**
     * @brief `at` returns a value specified by `index`. Note that this function does not protect the boundary.
     * @param index to return value.
     * @return double 
     */
    inline double at(int index) const {
        return data[index];
    }

    /**
     * @brief Set the `value` at the `index`. Note that this function does not protect the boundary.
     * 
     * @param index to set value at
     * @param value to set
     */
    inline void setAt(int index, double value) {
        data[index] = value;
    }
};

} }

#endif //UZUME_VOCODER_ROOT_WAVEFORM_HPP
