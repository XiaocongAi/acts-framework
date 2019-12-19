// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include <Acts/Utilities/Units.hpp>

#include "ACTFW/Digitization/HitSmearing.hpp"
#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/Fitting/FittingAlgorithm.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ACTFW/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Root/RootTrajectoryWriter.hpp"
#include "ACTFW/TruthTracking/ParticleSmearing.hpp"
#include "ACTFW/TruthTracking/TruthTrackFinder.hpp"
#include "ACTFW/Utilities/Options.hpp"
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
  for (auto cdr : geometry.second) { sequencer.addContextDecorator(cdr); }
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
  // TODO read truth hits
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

  // TODO pre-select particles

  // The fitter needs the measurements (proto tracks) and initial track states
  // (proto states). The elements in both collections must match and must be
  // created from the same input particles. Create truth tracks
  TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles       = particleReaderCfg.outputParticles;
  trackFinderCfg.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  trackFinderCfg.outputProtoTracks    = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));
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

  // setup the fitter
  FittingAlgorithm::Config fitCfg;
  fitCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  fitCfg.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  fitCfg.inputInitialTrackParameters
      = particleSmearingCfg.outputTrackParameters;
  fitCfg.outputTrajectories = "trajectories";
  fitCfg.fit                = FittingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(std::make_shared<FittingAlgorithm>(fitCfg, logLevel));

  // write tracks from fitting
  Root::RootTrajectoryWriter::Config trackWriterCfg;
  trackWriterCfg.inputParticles    = particleReaderCfg.outputParticles;
  trackWriterCfg.inputTrajectories = fitCfg.outputTrajectories;
  trackWriterCfg.outputDir         = outputDir;
  trackWriterCfg.outputFilename    = "tracks.root";
  trackWriterCfg.outputTreename    = "tracks";
  sequencer.addWriter(
      std::make_shared<Root::RootTrajectoryWriter>(trackWriterCfg, logLevel));
  // write reconstruction performance data
  TrackFinderPerformanceWriter::Config perFindCfg;
  perFindCfg.inputParticles       = particleReaderCfg.outputParticles;
  perFindCfg.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  perFindCfg.inputProtoTracks     = trackFinderCfg.outputProtoTracks;
  perFindCfg.outputDir            = outputDir;
  TrackFitterPerformanceWriter::Config perFitCfg;
  perFitCfg.inputParticles    = particleReaderCfg.outputParticles;
  perFitCfg.inputTrajectories = fitCfg.outputTrajectories;
  perFitCfg.outputDir         = outputDir;
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(perFindCfg, logLevel));
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perFitCfg, logLevel));

  return sequencer.run();
}
