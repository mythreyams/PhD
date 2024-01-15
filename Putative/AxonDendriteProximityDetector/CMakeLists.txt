# This project is designed to be built outside the Insight source tree.
PROJECT(AxonDendriteProximityFinder)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

LINK_LIBRARIES(ITKBasicFilters ITKIO )

#CMake Optimizing
SET(CMAKE_BUILD_TYPE RELEASE )
SET(CMAKE_CXX_FLAGS_DEBUG "-dDEBUG")
SET(CMAKE_C_FLAGS_DEBUG "-dDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(src)

SET(AxonDendriteProximityFinderSources src/amiraReader.h src/amiraReader.cpp src/proximity.h src/proximity.cpp src/axon_dendrite_proximity_finder.h src/axon_dendrite_proximity_finder.cpp src/axon_dendrite_proximity_finder_client.cpp src/basics.h src/basics.cpp src/bounding_box.h src/edge_wrapper.h src/typedefs.h src/typedefs_two.h src/utility.h src/utility.cpp src/brick_image_reader.h src/brick_image_reader.cpp src/pre_decon.cxx src/pre_decon.h)




ADD_EXECUTABLE( AxonDendriteProximityFinder ${AxonDendriteProximityFinderSources} )

TARGET_LINK_LIBRARIES(AxonDendriteProximityFinder ITKCommon ITKIO)
SET_TARGET_PROPERTIES(AxonDendriteProximityFinder PROPERTIES LINK_FLAGS "-openmp")
