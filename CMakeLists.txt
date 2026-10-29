cmake_minimum_required(VERSION 3.14)

project(DataPreprocessing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wl,--no-as-needed")

find_package(
  eudaq
  2.5.2
  COMPONENTS
  stave
  detector_event
  REQUIRED
)
include(${eudaq_DIR}/eudaqStaveConfig.cmake)

set(ISFEB2025 -1)
set(ISMARCH2026 -1)

string(FIND ${eudaq_DIR} "feb_2025" ISFEB2025)
string(FIND ${eudaq_DIR} "march_2026" ISMARCH2026)
if (${ISFEB2025} EQUAL -1 AND ${ISMARCH2026} EQUAL -1)
  message(FATAL_ERROR "eudaq_DIR is set to unknown location: " ${eudaq_DIR})
endif()

find_package(
  Eigen3 REQUIRED
)

find_package(
    ROOT 
    COMPONENTS 
    Core
    Hist
    HistPainter
    RIO 
    Tree
    Physics
    REQUIRED
)

find_package(
  OpenMP 
  COMPONENTS
  CXX
  REQUIRED
)

# Dictionary for the data
add_library(DetEventClusteredDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Types
)
ROOT_GENERATE_DICTIONARY(
    G__DetEventClusteredDict 
    ./include/Types/DetectorEvent.hpp
    MODULE 
    DetEventClusteredDict 
    LINKDEF 
    ./include/Types/LinkDef.hpp
)
target_link_libraries(
    DetEventClusteredDict
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
    eudaq::stave
    eudaq::detector_event
)

# Dictionary for the sim data
add_library(SimEventDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/detail
)
ROOT_GENERATE_DICTIONARY(
    G__SimEventDict
    ./include/detail/ROOTDictDefs.h
    MODULE 
    SimEventDict 
    LINKDEF 
    ./include/detail/LinkDef.h
)
target_link_libraries(
    SimEventDict 
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
)

if (NOT ${ISFEB2025} EQUAL -1)
  add_executable(
      dataPreprocessing 
      mainMembers.cpp
  )
elseif(NOT ${ISMARCH2026} EQUAL -1)
  add_executable(
      dataPreprocessing 
      mainMap.cpp
  )
endif()

add_executable(
    computeROIs 
    computeROIs.cpp
)

add_executable(
    simPreprocessing 
    mainSim.cpp
)

target_include_directories(
    dataPreprocessing
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/module/inc
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/hardware/inc
 )

target_include_directories(
    computeROIs
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/module/inc
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/hardware/inc
 )

target_include_directories(
    simPreprocessing
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/module/inc
    ${EUDAQ_INCLUDE_DIRS}/../user/stave/hardware/inc
 )

add_subdirectory(src/Analysis)
add_subdirectory(src/ROI)
add_subdirectory(src/detail)

target_link_libraries(
    dataPreprocessing
    PUBLIC
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    ROOT::Core
    ROOT::Hist
    ROOT::HistPainter
    ROOT::RIO
    ROOT::Tree
    ROOT::Physics
    eudaq::stave
    eudaq::detector_event
    DetEventClusteredDict 
)

target_link_libraries(
  computeROIs
  PUBLIC
  OpenMP::OpenMP_CXX
  Eigen3::Eigen
  ROOT::Core
  ROOT::Hist
  ROOT::HistPainter
  ROOT::RIO
  ROOT::Tree
  ROOT::Physics
  eudaq::stave
  eudaq::detector_event
  DetEventClusteredDict 
)

target_link_libraries(
  simPreprocessing
  PUBLIC
  ROOT::Core
  ROOT::Hist
  ROOT::HistPainter
  ROOT::RIO
  ROOT::Tree
  ROOT::Physics
  SimEventDict
)
