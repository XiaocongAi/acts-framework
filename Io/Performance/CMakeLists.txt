add_library(
  ActsFrameworkIoPerformance SHARED
  ACTFW/Io/Performance/TrackFinderPerformanceWriter.cpp
  ACTFW/Io/Performance/TrackFitterPerformanceWriter.cpp
  ACTFW/Io/Performance/CKFPerformanceWriter.cpp)
target_include_directories(
  ActsFrameworkIoPerformance
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsFrameworkIoPerformance
  PUBLIC ACTFramework
  PRIVATE ActsCore ROOT::Core ROOT::Tree)
