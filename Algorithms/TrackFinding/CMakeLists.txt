add_library(
  ACTFWTrackFinding SHARED
    src/TrackFindingAlgorithm.cpp
    src/TrackFindingAlgorithmTrackFinderFunction.cpp)
target_include_directories(
  ACTFWTrackFinding
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ACTFWTrackFinding
  PUBLIC ActsCore ACTFramework ACTFWBFieldPlugin Boost::program_options)

install(
  TARGETS ACTFWTrackFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
