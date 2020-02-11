// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/FakeRatePlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "ACTFW/Validation/TrackSummaryPlotTool.hpp"

class TFile;
class TTree;

namespace FW {

/// Write out the residual and pull of track parameters and efficiency.
///
/// Efficiency here is the fraction of smoothed tracks compared to all tracks.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class CKFPerformanceWriter final : public WriterT<CKFTrajectoryContainer>
{
public:
  struct Config
  {
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input (fitted) trajectories collection.
    std::string inputTrajectories;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_ckf.root";
    /// Plot tool configurations.
    ResPlotTool::Config          resPlotToolConfig;
    EffPlotTool::Config          effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    FakeRatePlotTool::Config     fakeRatePlotToolConfig;
  };

  /// Construct from configuration and log level.
  CKFPerformanceWriter(const Config& cfg, Acts::Logging::Level lvl);
  ~CKFPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode
  endRun() final override;

private:
  ProcessCode
  writeT(const AlgorithmContext&       ctx,
         const CKFTrajectoryContainer& trajectories) final override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile*     m_outputFile{nullptr};
  /// Plot tool for residuals and pulls.
  ResPlotTool               m_resPlotTool;
  ResPlotTool::ResPlotCache m_resPlotCache;
  /// Plot tool for efficiency
  EffPlotTool               m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache{};
  /// Plot tool for fake rate
  FakeRatePlotTool                    m_fakeRatePlotTool;
  FakeRatePlotTool::FakeRatePlotCache m_fakeRatePlotCache{};
  /// Plot tool for track summary
  TrackSummaryPlotTool                        m_trackSummaryPlotTool;
  TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache{};
};

}  // namespace FW
