# Abstract

Uzume is a simple wrapper library of WORLD's synthesis.
(WORLD is a high-quality speech analysis, manipulation and synthesis system by mmorise.
 See more: https://github.com/mmorise/World)

# License

Under MIT License. See LICENSE file.

# Usage

Write dependencies in your CMakeLists.txt:
```
## Dependencies
ExternalProject_add(uzume_vocoder
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/uzume_vocoder
        GIT_REPOSITORY https://github.com/haruneko/uzume_vocoder
        GIT_TAG master
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
        )
add_dependencies(your_app uzume_vocoder)
link_directories(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_include_directories(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib/libuzume_vocoder.a)
```

and then use uzume_vocoder in main.cpp:

```
#include <stdio.h>
#include <data/Waveform.hpp>
#include <spectrogram/StretchedPartialSpectrogram.hpp>
#include <spectrogram/WaveformSpectrogram.hpp>
#include <world/SynthesizeWaveformWithWORLD.hpp>

using namespace uzume::vocoder;
using namespace uzume::vocoder::world;

class SimpleDoubledTimeAxisMap : public uzume::vocoder::TimeAxisMap {
public:
    SimpleDoubledTimeAxisMap() = delete;
    SimpleDoubledTimeAxisMap(double msLength) : _msLength(msLength) { }
    ~SimpleDoubledTimeAxisMap() = default;
    double at(double ms) const { return ms / 2.0; }
    double msLength() const { return _msLength; }
private:
    const double _msLength;
};

int main(int argc, char *argv[]) {
    if(argc < 2 || 3 < argc) {
        printf("usage: uzume_vocoder_sample (in filepath) (out filepath)<optional>");
        exit(-1);
    }
    const char *inPath = argv[1];
    const char *outPath = argc == 3 ? argv[2] : "output.wav";

    // simply analyze spectrogram from waveform.
    auto *in = Waveform::read(inPath);
    auto *inSpec = new WaveformSpectrogram(in);

    // prepare time axis map, output spectrogram and output waveform.
    auto *tam = new SimpleDoubledTimeAxisMap(inSpec->msLength() * 2.0);
    auto *outSpec = new StretchedPartialSpectrogram(inSpec, tam);
    auto *out = new Waveform(in->length * 2, in->samplingFrequency);

    SynthesizeWaveformWithWORLD synthesize;

    // synthesize spectrogram in `outSpec` into output waveform in `out`.
    synthesize(out, outSpec);

    // save waveform.
    out->save(outPath);

    return 0;
}
```

See more detail in sample directory.
