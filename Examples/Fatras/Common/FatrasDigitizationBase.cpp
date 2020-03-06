// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasDigitizationBase.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Digitization/DigitizationAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Io/Root/RootPlanarClusterWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"

void
FW::setupDigitization(
    FW::Options::Variables&                       vars,
    FW::Sequencer&                                sequencer,
    std::shared_ptr<const FW::RandomNumbers>      randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry)
{
  // Read the standard options
  auto logLevel = FW::Options::readLogLevel(vars);

  // Set the module stepper
  auto pmStepper = std::make_shared<Acts::PlanarModuleStepper>(
      Acts::getDefaultLogger("PlanarModuleStepper", logLevel));

  // Read the digitization configuration
  FW::DigitizationAlgorithm::Config digiConfig;
  digiConfig.inputSimulatedHits  = "hits";
  digiConfig.outputClusters      = "clusters";
  digiConfig.planarModuleStepper = pmStepper;
  digiConfig.randomNumbers       = randomNumbers;
  digiConfig.trackingGeometry    = trackingGeometry;

  // Create the algorithm and add it to the sequencer
  sequencer.addAlgorithm(
      std::make_shared<FW::DigitizationAlgorithm>(digiConfig, logLevel));

  // Output directory
  std::string outputDir = vars["output-dir"].template as<std::string>();

  // Write digitisation output as Csv files
  if (vars["output-csv"].template as<bool>()) {
    // clusters as root
    FW::CsvPlanarClusterWriter::Config clusterWriterCsvConfig;
    clusterWriterCsvConfig.inputClusters      = digiConfig.outputClusters;
    clusterWriterCsvConfig.inputSimulatedHits = digiConfig.inputSimulatedHits;
    clusterWriterCsvConfig.outputDir          = outputDir;
    auto clusteWriterCsv = std::make_shared<FW::CsvPlanarClusterWriter>(
        clusterWriterCsvConfig, logLevel);
    // Add to the sequencer
    sequencer.addWriter(clusteWriterCsv);
  }

  // Write digitsation output as ROOT files
  if (vars["output-root"].template as<bool>()) {
    // clusters as root
    FW::RootPlanarClusterWriter::Config clusterWriterRootConfig;
    clusterWriterRootConfig.inputClusters      = digiConfig.outputClusters;
    clusterWriterRootConfig.inputSimulatedHits = digiConfig.inputSimulatedHits;
    clusterWriterRootConfig.filePath
        = FW::joinPaths(outputDir, digiConfig.outputClusters + ".root");
    clusterWriterRootConfig.treeName = digiConfig.outputClusters;
    auto clusteWriterRoot = std::make_shared<FW::RootPlanarClusterWriter>(
        clusterWriterRootConfig, logLevel);
    // Add to the sequencer
    sequencer.addWriter(clusteWriterRoot);
  }
}
