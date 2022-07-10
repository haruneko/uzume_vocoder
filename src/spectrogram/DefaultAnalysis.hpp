// Copyright 2022 Hal@shurabaP.  All rights reserved.
// Use of this source code is governed by a MIT style
// license that can be found in the LICENSE file.
#ifndef UZUME_VOCODER_DEFAULT_ANALYSIS_HPP
#define UZUME_VOCODER_DEFAULT_ANALYSIS_HPP

#include <functional>
#include "../AnalyzeAperiodicity.hpp"
#include "../AnalyzePeriodicity.hpp"
#include "../EstimateF0.hpp"

namespace uzume {
namespace vocoder {
    namespace DefaultAnalysis {
        extern const std::function<AnalyzeAperiodicity *(unsigned int)> AperiodicAnalysisFactory;
        extern const std::function<AnalyzePeriodicity *(unsigned int)> PeriodicAnalysisFactory;
        extern const std::function<EstimateF0 *(double)> F0EstimationFactory;
    }
}
}


#endif // UZUME_VOCODER_DEFAULT_ANALYSIS_HPP
