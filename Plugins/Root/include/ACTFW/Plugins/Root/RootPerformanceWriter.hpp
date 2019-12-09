// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

class TFile;
class TTree;

namespace FW {
namespace Root {

  using Identifier  = Data::SimSourceLink;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;
  using TrackState       = Acts::TrackState<Identifier, Acts::BoundParameters>;
  using TrajectoryVector = std::vector<
      std::pair<size_t, Acts::MultiTrajectory<Data::SimSourceLink>>>;

  /// @class RootPerformanceWriter
  ///
  /// Write out the residual and pull of track parameters and efficiency, i.e.
  /// fraction of smoothed track
  ///  into output file
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  class RootPerformanceWriter final : public WriterT<TrajectoryVector>
  {
  public:
    struct Config
    {
      std::string trackCollection;           ///< trajectory collection to write
      std::string simulatedEventCollection;  ///< truth particle collection
      std::string filePath;                  ///< path of the output file
      std::string fileMode = "RECREATE";     ///< file access mode
      TFile*      rootFile = nullptr;        ///< common root file
      ResPlotTool::Config resPlotToolConfig;  ///< ResPlotTool config
      EffPlotTool::Config effPlotToolConfig;  ///< EffPlotTool config
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    RootPerformanceWriter(const Config&        cfg,
                          Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootPerformanceWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// @brief Write method called by the base class
    /// @param [in] ctx is the algorithm context for event information
    /// @param [in] trajectories are what to be written out
    ProcessCode
    writeT(const AlgorithmContext& ctx,
           const TrajectoryVector& trajectories) final override;

  private:
    Config     m_cfg;         ///< The config class
    std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
    TFile*     m_outputFile{nullptr};  ///< The output file
    int        m_eventNr{0};           ///< The event number
    ResPlotTool::ResPlotCache
        m_resPlotCache;  ///< The cache object for residual/pull plots
    EffPlotTool::EffPlotCache
                 m_effPlotCache;  ///< The cache object for efficiency plots
    ResPlotTool* m_resPlotTool{nullptr};  ///< The plot tool for residual&pull
    EffPlotTool* m_effPlotTool{nullptr};  ///< The plot tool for efficiency
  };

}  // namespace Root
}  // namespace FW
