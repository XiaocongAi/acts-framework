// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <string>
#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Common/CollectionsOptions.hpp"
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/ParticleGun/ParticleGun.hpp"
#include "ACTFW/ParticleGun/ParticleGunOptions.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include "ACTFW/Random/RandomNumbersOptions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"

/// The main executable
///
/// Creates an instance of the particle gun algorithm
/// run particle generation and writes out

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // Add the standard options
  FW::Options::addCommonOptions<po::options_description>(desc, 1, 2);
  // Add the collections options/names
  FW::Options::addCollectionsOptions<po::options_description>(desc);
  // Add the random number options
  FW::Options::addRandomNumbersOptions<po::options_description>(desc);
  // Add the particle gun options
  FW::Options::addParticleGunOptions<po::options_description>(desc);
  // Add specific options for this example
  desc.add_options()("output-dir",
                     po::value<std::string>()->default_value(""),
                     "Output directory location.")(
      "output-root-file",
      po::value<std::string>()->default_value(""),
      "If the string is not empty: write a '.root' output file (full run).")(
      "output-csv-file",
      po::value<std::string>()->default_value(""),
      "If the string is not empty: Write a '.csv' output files (per event)");
  // map to store the given program options
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // Print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // Now read the common options
  auto logLevel = FW::Options::readLogLevel<po::variables_map>(vm);
  auto nEvents  = FW::Options::readNumberOfEvents<po::variables_map>(vm);

  // The evgen collection name
  std::string evgenCollection
      = FW::Options::readEvgenCollectionName<po::variables_map>(vm);

  // Create the random number engine
  auto randomNumberSvcCfg
      = FW::Options::readRandomNumbersConfig<po::variables_map>(vm);
  auto randomNumberSvc
      = std::make_shared<FW::RandomNumbersSvc>(randomNumberSvcCfg);

  // Create the barcode service
  FW::BarcodeSvc::Config barcodeSvcCfg;
  auto                   barcodeSvc = std::make_shared<FW::BarcodeSvc>(
      barcodeSvcCfg, Acts::getDefaultLogger("BarcodeSvc", logLevel));

  // Now read the particle gun configs
  auto particleGunCfg
      = FW::Options::readParticleGunConfig<po::variables_map>(vm);
  particleGunCfg.barcodeSvc      = barcodeSvc;
  particleGunCfg.randomNumberSvc = randomNumberSvc;
  particleGunCfg.evgenCollection = evgenCollection;
  auto particleGun               = std::make_shared<FW::ParticleGun>(
      particleGunCfg, Acts::getDefaultLogger("ParticleGun", logLevel));

  // Output directory
  std::string outputDir = vm["output-dir"].as<std::string>();

  // Write particles as CSV files
  std::shared_ptr<FW::Csv::CsvParticleWriter> pWriterCsv = nullptr;
  std::string csvFileName = vm["output-csv-file"].as<std::string>();
  if (!csvFileName.empty()) {
    if (csvFileName.find(".csv") == std::string::npos) {
      csvFileName += ".csv";
    }
    FW::Csv::CsvParticleWriter::Config pWriterCsvConfig;
    pWriterCsvConfig.collection     = evgenCollection;
    pWriterCsvConfig.outputDir      = outputDir;
    pWriterCsvConfig.outputFileName = csvFileName;
    pWriterCsv = std::make_shared<FW::Csv::CsvParticleWriter>(pWriterCsvConfig);
  }

  // Write particles as ROOT files
  std::shared_ptr<FW::Root::RootParticleWriter> pWriterRoot = nullptr;
  std::string rootFileName = vm["output-root-file"].as<std::string>();
  if (!rootFileName.empty()) {
    if (rootFileName.find(".root") == std::string::npos) {
      rootFileName += ".root";
    }
    // Write particles as ROOT TTree
    FW::Root::RootParticleWriter::Config pWriterRootConfig;
    pWriterRootConfig.collection = evgenCollection;
    pWriterRootConfig.filePath   = rootFileName;
    pWriterRootConfig.barcodeSvc = barcodeSvc;
    pWriterRoot
        = std::make_shared<FW::Root::RootParticleWriter>(pWriterRootConfig);
  }
  // create the config object for the sequencer
  FW::Sequencer::Config seqConfig;
  // now create the sequencer & add the relevant
  FW::Sequencer sequencer(seqConfig);
  sequencer.addServices({randomNumberSvc});
  sequencer.addReaders({particleGun});
  sequencer.appendEventAlgorithms({});
  if (pWriterRoot) sequencer.addWriters({pWriterRoot});
  if (pWriterCsv) sequencer.addWriters({pWriterCsv});
  sequencer.run(nEvents);
}
