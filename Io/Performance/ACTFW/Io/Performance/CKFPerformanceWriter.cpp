// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <numeric>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Io/Performance/CKFPerformanceWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"

FW::CKFPerformanceWriter::CKFPerformanceWriter(
    const FW::CKFPerformanceWriter::Config& cfg,
    Acts::Logging::Level                    lvl)
  : WriterT(cfg.inputTrajectories, "CKFPerformanceWriter", lvl)
  , m_cfg(cfg)
  , m_outputFile(
        TFile::Open(joinPaths(cfg.outputDir, cfg.outputFilename).c_str(),
                    "RECREATE"))
  , m_resPlotTool(m_cfg.resPlotToolConfig, lvl)
  , m_effPlotTool(m_cfg.effPlotToolConfig, lvl)
  , m_fakeRatePlotTool(m_cfg.fakeRatePlotToolConfig, lvl)
  , m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl)
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

  // initialize the plot tools
  m_resPlotTool.book(m_resPlotCache);
  m_effPlotTool.book(m_effPlotCache);
  m_fakeRatePlotTool.book(m_fakeRatePlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

FW::CKFPerformanceWriter::~CKFPerformanceWriter()
{
  m_resPlotTool.clear(m_resPlotCache);
  m_effPlotTool.clear(m_effPlotCache);
  m_fakeRatePlotTool.clear(m_fakeRatePlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);
  if (m_outputFile) { m_outputFile->Close(); }
}

FW::ProcessCode
FW::CKFPerformanceWriter::endRun()
{
  // fill residual and pull details into additional hists
  m_resPlotTool.refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool.write(m_resPlotCache);
    m_effPlotTool.write(m_effPlotCache);
    m_fakeRatePlotTool.write(m_fakeRatePlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::CKFPerformanceWriter::writeT(const AlgorithmContext&       ctx,
                                 const CKFTrajectoryContainer& trajectories)
{
  // Read truth particles from input collection
  const auto& particles
      = ctx.eventStore.get<SimParticles>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Containers holding the (non)reconstructed truth particles w/ count
  std::map<Barcode, size_t> matched{};
  std::map<Barcode, size_t> unmatched{};

  // Loop over all trajectories
  for (const CKFTrack& traj : trajectories) {
    if (!traj.hasTrajectory()) { continue; }

    const auto& [tips, track] = traj.trajectory();
    // -> std::pair<std::vector<size_t>, Acts::MultiTrajectory

    for (const size_t& tip : tips) {

      std::vector<ParticleHitCount> particleHitCount
          = traj.identifyMajorityParticle(tip);
      if (particleHitCount.empty()) { continue; }

      Barcode barcode_maj_truth_part = particleHitCount.front().particleId;
      auto    truth_part             = particles.find(barcode_maj_truth_part);
      // boost::container::vec_iterator<FW::Data::SimParticle *, true>

      // if (truth_part != particles.end()) {...

      size_t nStates       = 0;
      size_t nMeasurements = 0;
      size_t nOutliers     = 0;
      size_t nHoles        = 0;

      track.visitBackwards(tip, [&](const auto& state) {
        nStates++;
        auto typeFlags = state.typeFlags();
        if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
          nMeasurements++;
        } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
          nOutliers++;
        } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
          nHoles++;
        }
      });

      bool  is_matched = false;
      auto* target     = &unmatched;
      if (particleHitCount.front().hitCount * 1. / nMeasurements > .8) {
        is_matched = true;
        target     = &matched;
      }

      // count how often the truth particle was reconstructed
      if (target->count(barcode_maj_truth_part) == 0) {
        target->emplace(barcode_maj_truth_part, 0);
      }
      (*target)[barcode_maj_truth_part] += 1;

      // fill track summary plot
      // n{States,Measurements,Outliers,Holes}_vs_{pT,eta}
      m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache,
                                  *truth_part,
                                  nStates,
                                  nMeasurements,
                                  nOutliers,
                                  nHoles);

      // fill fakerate_vs_{pT,eta,phi} plots
      // fake rate = N_{reco, matched} / N_reco
      m_fakeRatePlotTool.fill(m_fakeRatePlotCache, *truth_part, !is_matched);

      // fill the residual plots
      if (traj.hasTrackParameters(tip)) {
        m_resPlotTool.fill(m_resPlotCache,
                           ctx.geoContext,
                           *truth_part,
                           traj.trackParameters(tip));
      }

    }  // <- loop over sub-trajectories in CKFTrack

  }  // <- loop over trajectories

  for (const auto& truth_particle : particles) {
    auto barcode = truth_particle.barcode();

    size_t fake_nTruthMatchedTracks = 0;
    size_t fake_nFakeTracks         = 0;

    auto it = matched.find(barcode);
    if (it != matched.end()) {  // check if reconstructed
      // efficiency = N_{truth matched, reco} / N_truth
      m_effPlotTool.fill(m_effPlotCache, truth_particle, true);
      fake_nTruthMatchedTracks = it->second;
    } else {
      m_effPlotTool.fill(m_effPlotCache, truth_particle, false);
    }

    it = unmatched.find(barcode);
    if (it != unmatched.end()) { fake_nFakeTracks = it->second; }

    m_fakeRatePlotTool.fill(m_fakeRatePlotCache,
                            truth_particle,
                            fake_nTruthMatchedTracks,
                            fake_nFakeTracks);
    // fills n{Reco,TruthMatched,Fake}Tracks and duplicationNum_vs_{pT,eta,phi}

  }  // <- loop over truth particles

  return ProcessCode::SUCCESS;
}
