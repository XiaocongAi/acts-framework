// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

#include "TFile.h"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

FW::Root::RootMaterialTrackWriter::RootMaterialTrackWriter(
    const FW::Root::RootMaterialTrackWriter::Config& cfg,
    Acts::Logging::Level                             level)
  : Base(cfg.collection, "RootMaterialTrackWriter", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  m_outputTree
      = new TTree(m_cfg.treeName.c_str(), "TTree from RootMaterialTrackWriter");
  if (m_outputTree == nullptr) throw std::bad_alloc();

  // Set the branches
  m_outputTree->Branch("v_x", &m_v_x);
  m_outputTree->Branch("v_y", &m_v_y);
  m_outputTree->Branch("v_z", &m_v_z);
  m_outputTree->Branch("v_px", &m_v_px);
  m_outputTree->Branch("v_py", &m_v_py);
  m_outputTree->Branch("v_pz", &m_v_pz);
  m_outputTree->Branch("v_phi", &m_v_phi);
  m_outputTree->Branch("v_eta", &m_v_eta);
  m_outputTree->Branch("t_X0", &m_tX0);
  m_outputTree->Branch("t_L0", &m_tL0);
  m_outputTree->Branch("mat_x", &m_step_x);
  m_outputTree->Branch("mat_y", &m_step_y);
  m_outputTree->Branch("mat_z", &m_step_z);
  m_outputTree->Branch("mat_step_length", &m_step_length);
  m_outputTree->Branch("mat_X0", &m_step_X0);
  m_outputTree->Branch("mat_L0", &m_step_L0);
  m_outputTree->Branch("mat_A", &m_step_A);
  m_outputTree->Branch("mat_Z", &m_step_Z);
  m_outputTree->Branch("mat_rho", &m_step_rho);

  if (m_cfg.prePostStep) {
    m_outputTree->Branch("mat_sx", &m_step_sx);
    m_outputTree->Branch("mat_sy", &m_step_sy);
    m_outputTree->Branch("mat_sz", &m_step_sz);
    m_outputTree->Branch("mat_ex", &m_step_ex);
    m_outputTree->Branch("mat_ey", &m_step_ey);
    m_outputTree->Branch("mat_ez", &m_step_ez);
  }
}

FW::Root::RootMaterialTrackWriter::~RootMaterialTrackWriter()
{
  m_outputFile->Close();
}

FW::ProcessCode
FW::Root::RootMaterialTrackWriter::endRun()
{
  // write the tree and close the file
  ACTS_INFO("Writing ROOT output File : " << m_cfg.filePath);
  m_outputFile->cd();
  m_outputTree->Write();
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootMaterialTrackWriter::writeT(
    const AlgorithmContext&                         ctx,
    const std::vector<Acts::RecordedMaterialTrack>& materialTracks)
{
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Loop over the material tracks and write them out
  for (auto& mtrack : materialTracks) {

    // Clearing the vector first
    m_step_sx.clear();
    m_step_sy.clear();
    m_step_sz.clear();
    m_step_x.clear();
    m_step_y.clear();
    m_step_z.clear();
    m_step_ex.clear();
    m_step_ey.clear();
    m_step_ez.clear();
    m_step_length.clear();
    m_step_X0.clear();
    m_step_L0.clear();
    m_step_A.clear();
    m_step_Z.clear();
    m_step_rho.clear();

    // Reserve the vector then
    size_t mints = mtrack.second.materialInteractions.size();
    m_step_sx.reserve(mints);
    m_step_sy.reserve(mints);
    m_step_sz.reserve(mints);
    m_step_x.reserve(mints);
    m_step_y.reserve(mints);
    m_step_ez.reserve(mints);
    m_step_ex.reserve(mints);
    m_step_ey.reserve(mints);
    m_step_ez.reserve(mints);
    m_step_length.reserve(mints);
    m_step_X0.reserve(mints);
    m_step_L0.reserve(mints);
    m_step_A.reserve(mints);
    m_step_Z.reserve(mints);
    m_step_rho.reserve(mints);

    // reset the global counter
    if (m_cfg.recalculateTotals) {
      m_tX0 = 0.;
      m_tL0 = 0.;
    } else {
      m_tX0 = mtrack.second.materialInX0;
      m_tL0 = mtrack.second.materialInL0;
    }

    // set the track information at vertex
    m_v_x   = mtrack.first.first.x();
    m_v_y   = mtrack.first.first.y();
    m_v_z   = mtrack.first.first.z();
    m_v_px  = mtrack.first.second.x();
    m_v_py  = mtrack.first.second.y();
    m_v_pz  = mtrack.first.second.z();
    m_v_phi = phi(mtrack.first.second);
    m_v_eta = eta(mtrack.first.second);

    // an now loop over the material
    for (auto& mint : mtrack.second.materialInteractions) {
      // The material step position information
      m_step_x.push_back(mint.position.x());
      m_step_y.push_back(mint.position.y());
      m_step_z.push_back(mint.position.z());

      if (m_cfg.prePostStep) {
        Acts::Vector3D prePos
            = mint.position - 0.5 * mint.pathCorrection * mint.direction;
        Acts::Vector3D posPos
            = mint.position + 0.5 * mint.pathCorrection * mint.direction;
        m_step_sx.push_back(prePos.x());
        m_step_sy.push_back(prePos.y());
        m_step_sz.push_back(prePos.z());
        m_step_ex.push_back(posPos.x());
        m_step_ey.push_back(posPos.y());
        m_step_ez.push_back(posPos.z());
      }
      // the material information
      const auto& mprops = mint.materialProperties;
      m_step_length.push_back(mprops.thickness());
      m_step_X0.push_back(mprops.averageX0());
      m_step_L0.push_back(mprops.averageL0());
      m_step_A.push_back(mprops.averageA());
      m_step_Z.push_back(mprops.averageZ());
      m_step_rho.push_back(mprops.averageRho());
      // re-calculate if defined to do so
      if (m_cfg.recalculateTotals) {
        m_tX0 += mprops.thicknessInX0();
        m_tL0 += mprops.thicknessInL0();
      }
    }
    // write to
    m_outputTree->Fill();
  }

  // return success
  return FW::ProcessCode::SUCCESS;
}
