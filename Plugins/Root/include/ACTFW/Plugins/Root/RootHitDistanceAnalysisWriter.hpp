#ifndef ACTFW_PLUGINS_ROOTHITDISTANCEANALYSISWRITER_H
#define ACTFW_PLUGINS_ROOTHITDISTANCEANALYSISWRITER_H

#include <mutex>
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Utilities/HitData.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Utilities/GeometryID.hpp"
#include "ACTS/Utilities/Logger.hpp"

class TFile;
class TTree;

namespace FW {

namespace Root {

  /// @class RootHitDistanceAnalysisWriter
  ///
  /// @brief Writes out the result of the hit distance analysis
  ///
  /// The RootHitDistanceAnalysisWriter writes out FW::AnalysisParameters for a
  /// given geometry object (layer/surface) + two additional parameters which
  /// can be used for e.g. position in r,z/eta,phi.
  /// This writer is intended to be used with the output from
  /// FW::HitDistanceAlgorithm.hpp .
  ///
  class RootHitDistanceAnalysisWriter
      : public FW::WriterT<std::map<Acts::GeometryID,
                                    std::tuple<FW::AnalysisParameters,
                                               FW::AnalysisParameters,
                                               double,
                                               double>>>
  {
  public:
    using Base = FW::WriterT<std::map<Acts::GeometryID,
                                      std::tuple<FW::AnalysisParameters,
                                                 FW::AnalysisParameters,
                                                 double,
                                                 double>>>;

    // @struct Config
    //
    // The nested config class
    struct Config
    {
    public:
      /// The input parameters of the distance analysis for component
      std::string hitAnalysis = "";
      /// The output file path
      std::string filePath = "hitAnalysis.root";
      /// The file access mode
      std::string fileMode = "RECREATE";
      /// The name of the output tree
      std::string treeName = "hitAnalysis";
    };

    /// Constructor
    /// @param cfg is the configuration class
    RootHitDistanceAnalysisWriter(const Config&        cfg,
                                  Acts::Logging::Level level
                                  = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootHitDistanceAnalysisWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// The protected writeT method, called by the WriterT base
    /// @param [in] ctx is the algorithm context for event consistency
    /// @param [in] params the analysis parameters to be written out
    ProcessCode
    writeT(
        const FW::AlgorithmContext& ctx,
        const std::map<Acts::GeometryID,
                       std::tuple<FW::AnalysisParameters,
                                  FW::AnalysisParameters,
                                  double,
                                  double>>& layerDistanceParams) final override;

    /// The config class
    Config m_cfg;
    /// Protect multi-threaded writes
    std::mutex m_writeMutex;
    /// The output file
    TFile* m_outputFile;
    /// The output tree
    TTree* m_outputTree;
    /// The GeometryID of the layer
    unsigned long long m_layerID;
    /// Optional parameter (e.g. for indicating position in r, eta)
    float m_par0;
    /// Optional parameter (e.g. for indicating position in z, phi)
    float m_par1;
    /// Mean of component 0
    float m_mean0;
    /// Minimum of component 0
    float m_min0;
    /// Maximum of component 0
    float m_max0;
    /// Mean of component 1
    float m_mean1;
    /// Minimum of component 1
    float m_min1;
    /// Maximum of component 1
    float m_max1;
  };

}  // namespace Root
}  // namespace FW

#endif  // ACTFW_PLUGINS_ROOTHITDISTANCEANALYSISWRITER_H