add_library(
  ACTFWFinding SHARED
  src/FindingAlgorithm.cpp
  src/FindingAlgorithmFinderFunction.cpp)
target_include_directories(
  ACTFWFinding
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ACTFWFinding
  PUBLIC ActsCore ACTFramework ACTFWBFieldPlugin Boost::program_options)

install(
  TARGETS ACTFWFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
