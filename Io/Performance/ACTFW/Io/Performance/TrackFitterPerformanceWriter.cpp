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

  for (const auto& traj : trajectories) {
    if (traj.empty()) { continue; }
    // retrieve the truth particle barcode for this track state
    // TODO use majority particle instead to support reconstructed, non-truth
    // trajectories
    auto truthHitAtFirstState = (*traj[0].measurement.uncalibrated).truthHit();
    auto barcode              = truthHitAtFirstState.particle.barcode();
    // find the truth particle for this trajectory
    Data::SimParticle truthParticle;
    auto              ip = particles.find(barcode);
    if (ip != particles.end()) {
      truthParticle = *ip;
    } else {
      ACTS_WARNING("Truth particle " << barcode << "not found.");
    }
    // fill the plots
    m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, traj);
    m_effPlotTool.fill(m_effPlotCache, traj, truthParticle);
  }

  return ProcessCode::SUCCESS;
}
