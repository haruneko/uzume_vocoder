// Copyright 2020 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_WAVEFORM_HPP
#define UZUME_VOCODER_WAVEFORM_HPP

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
    double msLength() const;

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
};

} }

#endif //UZUME_VOCODER_ROOT_WAVEFORM_HPP
