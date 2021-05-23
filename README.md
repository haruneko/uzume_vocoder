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
#include <algorithm>

#include "SynthesizeImpulseResponseWithWORLD.hpp"
#include "SynthesizeSegmentWithWORLD.hpp"
#include "Waveform.hpp"
#include "WaveformSpectrogram.hpp"

using namespace uzume::vocoder;

// This is a sample to use vocoder directory.
int main() {
    const char *inputPath = "/path/to/input.wav";
    const char *outputPath = "/path/to/output.wav";
    Waveform *input = Waveform::read(inputPath);
    Waveform *output = new Waveform(input->length, input->samplingFrequency);
    WaveformSpectrogram spectrogram(input);

    SynthesizeImpulseResponseWithWORLD irs(spectrogram.fftSize(), input->samplingFrequency);
    SynthesizeSegmentWithWORLD synthesize(&irs);

    for(unsigned int i = 0; i < output->length; i++) {
        output->data[i] = 0.0;
    }

    SegmentSignal s(output->data, /* indexMin = */ 0, /* indexMax = */ output->length, output->samplingFrequency);
    SegmentParameters p(&spectrogram, /* startPhase = */ 0.0, /* startFractionalTimeShift = */ 0.0);

    synthesize(&s, &p);

    output->save(outputPath);
}
```
