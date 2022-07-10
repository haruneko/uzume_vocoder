// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#include "DefaultAnalysis.hpp"

#include "../world/AnalyzeAperiodicityWithD4C.hpp"
#include "../world/AnalyzePeriodicityWithCheapTrick.hpp"
#include "../world/EstimateF0WithDIO.hpp"

using namespace uzume::vocoder;
using namespace uzume::vocoder::world;

const std::function<AnalyzeAperiodicity *(unsigned int)> DefaultAnalysis::AperiodicAnalysisFactory = [](unsigned int samplingFrequency) {
    return new AnalyzeAperiodicityWithD4C(samplingFrequency);
};

const std::function<AnalyzePeriodicity *(unsigned int)> DefaultAnalysis::PeriodicAnalysisFactory = [](unsigned int samplingFrequency) {
    return new AnalyzePeriodicityWithCheapTrick(samplingFrequency);
};

const std::function<EstimateF0 *(double)> DefaultAnalysis::F0EstimationFactory = [](double msFramePeriod) {
    return new EstimateF0WithDIO(msFramePeriod);
};
