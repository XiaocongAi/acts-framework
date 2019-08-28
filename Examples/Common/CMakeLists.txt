#file(GLOB_RECURSE src_files "src")
set(src_files src/CommonGeometry.cpp src/CommonOptions.cpp src/GeometryExampleBase.cpp src/MaterialMappingBase.cpp src/MaterialValidationBase.cpp src/PropagationExampleBase.cpp)

# append Fatras example is Pythia8 is avaialbe
if(USE_PYTHIA8)
  list(APPEND src_files src/FatrasExampleBase.cpp)
endif()

add_library(ACTFWExamplesCommon SHARED ${src_files})
target_include_directories(ACTFWExamplesCommon PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWExamplesCommon PUBLIC ActsCore FatrasCore ACTFramework)
target_link_libraries(ACTFWExamplesCommon PUBLIC ACTFWObjPlugin 
        ACTFWCsvPlugin ACTFWJsonPlugin ACTFWRootPlugin ACTFWDetectorsCommon
        ACTFWBFieldPlugin ACTFWDigitization ACTFWFatras ActsFrameworkGenerators ACTFWMaterialMapping)

if(USE_PYTHIA8)
  target_link_libraries(ACTFWExamplesCommon PUBLIC ActsFrameworkPythia8)
endif()

install(TARGETS ACTFWExamplesCommon LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})