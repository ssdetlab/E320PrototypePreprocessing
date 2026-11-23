if (NOT ${ISFEB2025} EQUAL -1)
  message("-- Found Feb 2025 eudaq implementation. Building members-based analysis")
  message(FATAL_ERROR "-- Feb 2025 implementation is deprecated")
elseif(NOT ${ISMARCH2026} EQUAL -1)
  message("-- Found March 2026 eudaq implementation. Building map-based analysis")
  target_sources(
      dataPreprocessing
      PRIVATE
      PreprocessorMap.cpp
      ClusteringAlgorithm.cpp
      AcceleratorStateFilter.cpp
  )
endif()

target_sources(
    simPreprocessing
    PRIVATE
    PreprocessorSim.cpp
    ClusteringAlgorithm.cpp
)
