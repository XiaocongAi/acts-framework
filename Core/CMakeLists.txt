add_library(ACTFramework SHARED
  src/Framework/BareAlgorithm.cpp
  src/Framework/BareService.cpp
  src/Framework/RandomNumbers.cpp
  src/Framework/Sequencer.cpp
  src/Framework/CalibrationDecorator.cpp
  src/Utilities/Paths.cpp
  src/Utilities/Helpers.cpp
  src/Validation/EffPlotTool.cpp
  src/Validation/FakeRatePlotTool.cpp
  src/Validation/TrackSummaryPlotTool.cpp
  src/Validation/ProtoTrackClassification.cpp
  src/Validation/ResPlotTool.cpp)
target_include_directories(
  ACTFramework
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(
  ACTFramework
  PUBLIC ActsCore Boost::boost ROOT::Core ROOT::Hist
  PRIVATE ${TBB_LIBRARIES} Boost::filesystem dfelibs)
target_compile_definitions(
  ACTFramework
  PRIVATE BOOST_FILESYSTEM_NO_DEPRECATED)
# set per-target c++17 requirement that will be propagated to linked targets
target_compile_features(ACTFramework PUBLIC cxx_std_17)

install(
  TARGETS ACTFramework
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
