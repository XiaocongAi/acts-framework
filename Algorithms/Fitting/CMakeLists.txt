add_library(
  ACTFWFitting SHARED
  src/FittingAlgorithm.cpp
  src/FittingAlgorithmFitterFunction.cpp)
target_include_directories(
  ACTFWFitting
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(
  ACTFWFitting
  PUBLIC ActsCore ACTFramework ACTFWBFieldPlugin Boost::program_options
  PRIVATE ${TBB_LIBRARIES})

install(
  TARGETS ACTFWFitting
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
