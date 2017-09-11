#ifndef ACTFW_BFIELD_BFIELDEXAMPLE_H
#define ACTFW_BFIELD_BFIELDEXAMPLE_H

#include <string>
#include <boost/program_options.hpp>
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/StandardOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/BField/RootInterpolatedBFieldWriter.hpp"

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file and writes out the
/// grid points and values of the map into root format. The Field can then be
/// displayed using the root script printBField.cpp

namespace po = boost::program_options;

int
main(int argc, char* argv[])
{
  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // add the standard options
  FW::Options::addStandardOptions<po::options_description>(desc,1,2);
  // add the bfield options
  FW::Options::addBFieldOptions<po::options_description>(desc);
  // add an output file
  desc.add_options()("out", 
          po::value<std::string>()->default_value("BFieldOut.root"),
          "Set this name for an output root file.");
            
  // map to store the given program options
  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // read the standard options
  // print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  // now read the standard options options
  auto standardOptions 
    = FW::Options::readStandardOptions<po::variables_map>(vm);
  auto nEvents = standardOptions.first;
  auto logLevel = standardOptions.second;
  // create BField service
  auto bField = FW::Options::readBField<po::variables_map>(vm);
  if (!bField.first) {
    std::cout << "Bfield could not be set up. Exiting." << std::endl;
    return -1;
  }
    
  // Create the InterpolatedBFieldWriter
  FW::BField::RootInterpolatedBFieldWriter::Config writerConfig;
  if (vm["rz"].as<bool>())
    writerConfig.gridType = FW::BField::GridType::rz;
  else
    writerConfig.gridType = FW::BField::GridType::xyz;
  writerConfig.treeName   = "bField";
  writerConfig.fileName   = vm["out"].as<std::string>();

  writerConfig.bField = bField.first;
  auto bFieldWriter
      = std::make_shared<FW::BField::RootInterpolatedBFieldWriter>(
          writerConfig);

  // create the config object for the sequencer
  FW::Sequencer::Config seqConfig;
  // now create the sequencer
  FW::Sequencer sequencer(seqConfig);
  sequencer.addServices({bFieldWriter});
  sequencer.run(1);
}

#endif  // ACTFW_BFIELD_BFIELDEXAMPLE_H
