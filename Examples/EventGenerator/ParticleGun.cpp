// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <memory>

#include <Acts/Utilities/Units.hpp>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Io/Root/RootParticleWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::units;
using namespace FW;

int
main(int argc, char* argv[])
{
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addParticleGunOptions(desc);
  Options::addOutputOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  Acts::Logging::Level logLevel = Options::readLogLevel(vm);

  // basic services
  auto rnd
      = std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // event generation w/ particle gun
  EventGenerator::Config evgenCfg = Options::readParticleGunOptions(vm);
  evgenCfg.output                 = "particles";
  evgenCfg.randomNumbers          = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgenCfg, logLevel));

  // different output modes
  std::string outputDir = vm["output-dir"].as<std::string>();
  if (vm["output-csv"].as<bool>()) {
    CsvParticleWriter::Config csvWriterCfg;
    csvWriterCfg.inputEvent = evgenCfg.output;
    csvWriterCfg.outputDir  = outputDir;
    csvWriterCfg.outputStem = "particles";
    sequencer.addWriter(
        std::make_shared<CsvParticleWriter>(csvWriterCfg, logLevel));
  }
  if (vm["output-root"].as<bool>()) {
    RootParticleWriter::Config rootWriterCfg;
    rootWriterCfg.collection = evgenCfg.output;
    rootWriterCfg.filePath   = joinPaths(outputDir, "particles.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(rootWriterCfg, logLevel));
  }

  return sequencer.run();
}
