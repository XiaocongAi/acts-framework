// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootPerformanceWriter.hpp"
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;

FW::Root::RootPerformanceWriter::RootPerformanceWriter(
    const FW::Root::RootPerformanceWriter::Config& cfg,
    Acts::Logging::Level                           level)
  : WriterT(cfg.trackCollection, "RootPerformanceWriter", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // Input track and truth collection name
  if (m_cfg.trackCollection.empty()) {
    throw std::invalid_argument("Missing input trajectory collection");
  } else if (m_cfg.simulatedEventCollection.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }

  // Initialize the residual and efficiency plots tool
  m_resPlotTool = new FW::ResPlotTool(m_cfg.resPlotToolConfig, level);
  m_effPlotTool = new FW::EffPlotTool(m_cfg.effPlotToolConfig, level);
  if (m_resPlotTool == nullptr) {
    throw std::bad_alloc();
  } else if (m_effPlotTool == nullptr) {
    throw std::bad_alloc();
  }

  m_resPlotTool->book(m_resPlotCache);
  m_effPlotTool->book(m_effPlotCache);
}

FW::Root::RootPerformanceWriter::~RootPerformanceWriter()
{
  m_resPlotTool->clear(m_resPlotCache);
  m_effPlotTool->clear(m_effPlotCache);
  delete m_resPlotTool;
  delete m_effPlotTool;
  if (m_outputFile) { m_outputFile->Close(); }
}

FW::ProcessCode
FW::Root::RootPerformanceWriter::endRun()
{
  // fill residual and pull details into additional hists
  m_resPlotTool->refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool->write(m_resPlotCache);
    m_effPlotTool->write(m_effPlotCache);
    ACTS_INFO("Write performance plots to '" << m_cfg.filePath << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootPerformanceWriter::writeT(const AlgorithmContext& ctx,
                                        const TrajectoryVector& trajectories)
{
  if (m_outputFile == nullptr) return ProcessCode::SUCCESS;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Read truth particles from input collection
  const auto& simulatedEvent = ctx.eventStore.get<std::vector<Data::SimVertex>>(
      m_cfg.simulatedEventCollection);
  ACTS_DEBUG("Read collection '" << m_cfg.simulatedEventCollection << "' with "
                                 << simulatedEvent.size() << " vertices");

  // Get the map of truth particle
  ACTS_DEBUG("Get the truth particles.");
  std::map<barcode_type, Data::SimParticle> particles;
  for (auto& vertex : simulatedEvent) {
    for (auto& particle : vertex.outgoing) {
      particles.insert(std::make_pair(particle.barcode(), particle));
    }
  }

  // Loop over the trajectories
  for (auto& [trackTip, traj] : trajectories) {
    // retrieve the truth particle barcode for this track state
    auto truthHitAtFirstState
        = traj.getTrackState(trackTip).uncalibrated().truthHit();

    traj.visitBackwards(trackTip, [&](const auto& state) {
      auto truthHitAtFirstState = state;
    });

    auto barcode = truthHitAtFirstState.particle.barcode();
    // find the truth Particle for this trajectory
    Data::SimParticle truthParticle;
    if (particles.find(barcode) != particles.end()) {
      ACTS_DEBUG("Find the truth particle with barcode = " << barcode);
      truthParticle = particles.find(barcode)->second;
    } else {
      ACTS_WARNING("Truth particle with barcode = " << barcode << "not found.");
    }

    // fill the plots
    m_resPlotTool->fill(m_resPlotCache, ctx.geoContext, traj, trackTip);
    m_effPlotTool->fill(m_effPlotCache, traj, trackTip, truthParticle);

  }  // all trajectories

  return ProcessCode::SUCCESS;
}
