add_executable(
    ActsRecTruthTracks
    RecTruthTracks.cpp)

target_link_libraries(
   ActsRecTruthTracks
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
   ActsRecTruthTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
