// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Performance/TrackFitterPerformanceWriter.hpp"

#include <stdexcept>

#include <Acts/Utilities/Helpers.hpp>
#include <TFile.h>
#include <TTree.h>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using Acts::VectorHelpers::eta;

FW::TrackFitterPerformanceWriter::TrackFitterPerformanceWriter(
    const FW::TrackFitterPerformanceWriter::Config& cfg,
    Acts::Logging::Level                            lvl)
  : WriterT(cfg.inputTrajectories, "TrackFitterPerformanceWriter", lvl)
  , m_cfg(cfg)
  , m_outputFile(
        TFile::Open(joinPaths(cfg.outputDir, cfg.outputFilename).c_str(),
                    "RECREATE"))
  , m_resPlotTool(m_cfg.resPlotToolConfig, lvl)
  , m_effPlotTool(m_cfg.effPlotToolConfig, lvl)
{
  // Input track and truth collection name
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path    = joinPaths(cfg.outputDir, cfg.outputFilename);
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the residual and efficiency plots tool
  m_resPlotTool.book(m_resPlotCache);
  m_effPlotTool.book(m_effPlotCache);
}

FW::TrackFitterPerformanceWriter::~TrackFitterPerformanceWriter()
{
  m_resPlotTool.clear(m_resPlotCache);
  m_effPlotTool.clear(m_effPlotCache);
  if (m_outputFile) { m_outputFile->Close(); }
}

FW::ProcessCode
FW::TrackFitterPerformanceWriter::endRun()
{
  // fill residual and pull details into additional hists
  m_resPlotTool.refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool.write(m_resPlotCache);
    m_effPlotTool.write(m_effPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::TrackFitterPerformanceWriter::writeT(
    const AlgorithmContext&    ctx,
    const TrajectoryContainer& trajectories)
{
  // read truth particles from input collection
  const auto& particles
      = ctx.eventStore.get<SimParticles>(m_cfg.inputParticles);

  // exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // store the reconstructed trajectories with truth info
  std::map<barcode_type, std::pair<size_t, Trajectory>> reconTrajectories;
  barcode_type                                          barcode{0};
  bool                                                  hasMeasurement = false;
  for (const auto& [trackTip, traj] : trajectories) {
    // retrieve the truth particle barcode for this trajectory
    traj.visitBackwards(trackTip, [&](const auto& state) {
      // only the track state with measurement has truth hit
      if (state.hasUncalibrated()) {
        auto truthHit = state.uncalibrated().truthHit();
        // retrieve the truth particle barcode for this track state
        barcode        = truthHit.particle.barcode();
        hasMeasurement = true;
        return false;  // abort the execution after we get the barcode
      }
      return true;  // continue the execution
    });

    // no further process if no measurements at all
    if (not hasMeasurement) {
      continue;
      ACTS_DEBUG("No measurements on this track!");
    }

    // record this trajectory with its truth info
    reconTrajectories.emplace(barcode, std::make_pair(trackTip, traj));

    // fill the residual plots
    m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, traj, trackTip);
  }

  // Fill the efficiency plots
  // @Todo: store the truth particle in the track collection even the fit fail
  // to simplify the filling of efficiency
  // The defintion of efficiency is W.R.T. all truth particles
  for (const auto& particle : particles) {
    auto barcode = particle.barcode();
    auto ip      = reconTrajectories.find(barcode);
    if (ip != reconTrajectories.end()) {
      const auto& [trackTip, traj] = ip->second;
      // when the trajectory is reconstructed
      m_effPlotTool.fill(m_effPlotCache, traj, trackTip, particle);
    } else {
      // when the trajectory is NOT reconstructed
      m_effPlotTool.fill(m_effPlotCache, particle);
    }
  }

  return ProcessCode::SUCCESS;
}
