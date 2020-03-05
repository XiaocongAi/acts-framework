add_executable(
    ACTFWRecTruthTracksExample
    RecTruthTracks.cpp)

target_link_libraries(
    ACTFWRecTruthTracksExample
    PRIVATE
      ACTFramework
      ACTFWFitting
      ACTFWGenericDetector
      ACTFWExamplesCommon
      ACTFWBFieldPlugin
      ActsFrameworkTruthTracking
      ActsFrameworkIoCsv
      ActsFrameworkIoPerformance)

install(
  TARGETS
    ACTFWRecTruthTracksExample
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
