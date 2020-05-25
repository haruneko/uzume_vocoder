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
ExternalProject_add(uzume_dsp
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/uzume_dsp
        GIT_REPOSITORY https://github.com/haruneko/uzume_dsp
        GIT_TAG master
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
        )
add_dependencies(your_app uzume_dsp)
link_directories(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_include_directories(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(your_app PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib/libuzume_dsp.a)
```

and then use uzume_dsp in main.cpp:

```
#include <algorithm>

#include "SynthesizeImpulseResponseWithWORLD.hpp"
#include "SynthesizePhraseWithWORLD.hpp"
#include "Waveform.hpp"
#include "WaveformSpectrogram.hpp"

using namespace uzume::dsp;

// This is a sample to use dsp directory.
int main() {
    const char *inputPath = "/path/to/input.wav";
    const char *outputPath = "/path/to/output.wav";
    Waveform *input = Waveform::read(inputPath);
    Waveform *output = new Waveform(input->length, input->samplingFrequency);
    WaveformSpectrogram spectrogram(input);

    SynthesizeImpulseResponseWithWORLD irs(spectrogram.fftSize(), input->samplingFrequency);
    SynthesizePhraseWithWORLD synthesize(&irs);

    for(unsigned int i = 0; i < output->length; i++) {
        output->data[i] = 0.0;
    }

    PhraseSignal s(output->data, /* indexMin = */ 0, /* indexMax = */ output->length, output->samplingFrequency);
    PhraseParameters p(&spectrogram, /* startPhase = */ 0.0, /* startFractionalTimeShift = */ 0.0);

    synthesize(&s, &p);

    output->save(outputPath);
}
```
