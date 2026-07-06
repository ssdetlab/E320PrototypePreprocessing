add_library(eudaq_preprocessor::detector_event SHARED IMPORTED)
set_target_properties(eudaq_preprocessor::detector_event PROPERTIES
  IMPORTED_LOCATION "${CMAKE_CURRENT_LIST_DIR}/../build/libDetEventClusteredDict.so"
  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/../include/Types"
)
