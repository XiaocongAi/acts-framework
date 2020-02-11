// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

// #include <ACTFW/Io/Root/RootTrajectoryWriter.hpp>
#include <Acts/Utilities/Units.hpp>

#include "ACTFW/Digitization/HitSmearing.hpp"
#include "ACTFW/Finding/FindingAlgorithm.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Io/Performance/CKFPerformanceWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TruthTracking/ParticleSmearing.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::UnitLiterals;
using namespace FW;

int
main(int argc, char* argv[])
{
  GenericDetector detector;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  detector.addOptions(desc);
  Options::addBFieldOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel  = Options::readLogLevel(vm);
  auto inputDir  = vm["input-dir"].as<std::string>();
  auto outputDir = vm["output-dir"].as<std::string>();
  auto rnd       = std::make_shared<FW::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  // ensure the output directory exists
  ensureWritableDirectory(outputDir);

  // Setup detector geometry
  auto geometry         = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (const auto& cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);

  // Read particles and clusters from CSV files
  auto particleReaderCfg            = Options::readCsvParticleReaderConfig(vm);
  particleReaderCfg.outputParticles = "truth_particles";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReaderCfg, logLevel));
  // Read clusters from CSV files
  auto clusterReaderCfg = Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = trackingGeometry;

  clusterReaderCfg.outputClusters        = "clusters";
  clusterReaderCfg.outputHitIds          = "hit_ids";
  clusterReaderCfg.outputHitParticlesMap = "truth_hit_particles_map";
  clusterReaderCfg.outputSimulatedHits   = "truth_hits";
  sequencer.addReader(
      std::make_shared<CsvPlanarClusterReader>(clusterReaderCfg, logLevel));

  // Create smeared measurements
  HitSmearing::Config hitSmearingCfg;
  hitSmearingCfg.inputSimulatedHits = clusterReaderCfg.outputSimulatedHits;
  hitSmearingCfg.outputSourceLinks  = "sourcelinks";
  hitSmearingCfg.sigmaLoc0          = 25_um;
  hitSmearingCfg.sigmaLoc1          = 100_um;
  hitSmearingCfg.randomNumbers      = rnd;
  sequencer.addAlgorithm(
      std::make_shared<HitSmearing>(hitSmearingCfg, logLevel));
  // Create smeared particles states
  ParticleSmearing::Config particleSmearingCfg;
  particleSmearingCfg.inputParticles        = particleReaderCfg.outputParticles;
  particleSmearingCfg.outputTrackParameters = "smearedparameters";
  particleSmearingCfg.randomNumbers         = rnd;
  // Gaussian sigmas to smear particle parameters
  particleSmearingCfg.sigmaD0    = 20_um;
  particleSmearingCfg.sigmaD0PtA = 30_um;
  particleSmearingCfg.sigmaD0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaZ0    = 20_um;
  particleSmearingCfg.sigmaZ0PtA = 30_um;
  particleSmearingCfg.sigmaZ0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaPhi   = 1_degree;
  particleSmearingCfg.sigmaTheta = 1_degree;
  particleSmearingCfg.sigmaPRel  = 0.01;
  particleSmearingCfg.sigmaT0    = 1_ns;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(particleSmearingCfg, logLevel));

  // setup finding algorithm
  FindingAlgorithm::Config findCfg;
  findCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  findCfg.inputInitialTrackParameters
      = particleSmearingCfg.outputTrackParameters;
  findCfg.outputTrajectories = "trajectories";
  findCfg.find               = FindingAlgorithm::makeFinderFunction(
      trackingGeometry, magneticField, logLevel);
  // findCfg.slsCfg = Acts::CKFSourceLinkSelector::Config(...);
  sequencer.addAlgorithm(std::make_shared<FindingAlgorithm>(findCfg, logLevel));

  // write tracks
  // TODO: RootTrajectoryWriter does not like CKFTrajectoryContainers
  // FW::RootTrajectoryWriter::Config trackWriterCfg;
  // trackWriterCfg.inputParticles    = particleReaderCfg.outputParticles;
  // trackWriterCfg.inputTrajectories = findCfg.outputTrajectories;
  // trackWriterCfg.outputDir         = outputDir;
  // trackWriterCfg.outputFilename    = "tracks.root";
  // trackWriterCfg.outputTreename    = "tracks";
  // sequencer.addWriter(
  //     std::make_shared<FW::RootTrajectoryWriter>(trackWriterCfg, logLevel));

  // write reconstruction performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles    = particleReaderCfg.outputParticles;
  perfWriterCfg.inputTrajectories = findCfg.outputTrajectories;
  perfWriterCfg.outputDir         = outputDir;
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  return sequencer.run();
}
