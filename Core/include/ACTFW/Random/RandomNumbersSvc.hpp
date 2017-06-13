//
//  RandomNumbersSvc.hpp
//  ACTFW
//
//  Created by Andreas Salzburger on 17/05/16.
//
//

#ifndef ACTFW_RANDOM_RANDOMNUMBERSSVC_H
#define ACTFW_RANDOM_RANDOMNUMBERSSVC_H 1

#include <array>
#include <memory>
#include <random>
#include <string>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace FW {

/// @class RandomNumbersSvc
///
/// This service provides Algorithm-local random number generators, allowing for
/// thread-safe, lock-free and reproducible random number generations in both
/// single-threaded and multi-threaded test framework runs.
///
typedef std::mt19937                     RandomEngine;  ///< Mersenne Twister
typedef std::normal_distribution<double> GaussDist;     ///< Normal Distribution
typedef std::uniform_real_distribution<double>
                                        UniformDist;  ///< Uniform Distribution
typedef std::gamma_distribution<double> GammaDist;    ///< Gamma Distribution
typedef std::poisson_distribution<int>  PoissonDist;  ///< Poisson Distribution

class RandomNumbersSvc : public IService
{
public:
  /// @class Config
  ///
  /// Nested Configuration class
  struct Config
  {
    /// default logger
    std::shared_ptr<const Acts::Logger> logger;
    /// service name
    std::string name;
    /// random seed
    unsigned int seed = 1234567890;
    /// configuration uniform
    std::array<double, 2> uniform_parameters = {{0, 1}};
    /// configuration gauss
    std::array<double, 2> gauss_parameters = {{0, 1}};
    /// configuration landau
    std::array<double, 2> landau_parameters = {{0, 1}};
    /// configuration gamma
    std::array<double, 2> gamma_parameters = {{0, 1}};
    /// configuration poisson
    int poisson_parameter = 40;

    Config(const std::string&   lname = "RandomNumbersSvc",
           Acts::Logging::Level lvl   = Acts::Logging::INFO)
      : logger{Acts::getDefaultLogger(lname, lvl)}, name{lname}
    {
    }
  };

  /// @class Generator
  ///
  /// A random number generator. The intended mode of operation is that each
  /// Algorithm::execute() invocation should get its own private instance of
  /// this class.
  ///
  class Generator
  {
  public:
    /// Initialize a generator
    ///
    /// @param cfg is the host's configuration
    /// @param seed is the seed to use for this generator
    Generator(const Config& cfg, unsigned int seed);
    
    /// draw random number from gauss distribution
    double
    drawGauss();
    /// draw random number from uniform distribution
    double
    drawUniform();
    /// draw random number from landau distribution
    double
    drawLandau();
    /// draw random number from gamma distribution
    double
    drawGamma();
    /// draw random number from poisson distribution
    double
    drawPoisson();

  private:
    const Config& m_cfg;      ///< link to host configuration
    RandomEngine  m_engine;   ///< random engine
    GaussDist     m_gauss;    ///< gauss distribution
    UniformDist   m_uniform;  ///< uniform distribution
    GammaDist     m_gamma;    ///< gamma distribution
    PoissonDist   m_poisson;  ///< poisson distribution
  };

  /// Constructor
  RandomNumbersSvc(const Config& cfg);

  // Framework initialize method
  FW::ProcessCode
  initialize() override final;

  /// Framework finalize mehtod
  FW::ProcessCode
  finalize() override final;

  /// Spawn an algorithm-local random number generator
  ///
  /// @param context is the AlgorithmContext of the host algorithm
  Generator
  spawnGenerator(const AlgorithmContext& context) const;

  /// Ask for the seed
  unsigned int
  seed() const
  {
    return m_cfg.seed;
  }

  /// Framework name() method
  const std::string&
  name() const override final;

private:
  Config    m_cfg;  ///< the configuration class

  /// Private access to the logging instance
  const Acts::Logger&
  logger() const
  {
    return *m_cfg.logger;
  }
};

inline const std::string&
RandomNumbersSvc::name() const
{
  return m_cfg.name;
}
}

#endif  // ACTFW_RANDOM_RANDOMNUMBERSSVC_H
