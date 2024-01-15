# This project is designed to be built outside the Insight source tree.
PROJECT(BoutonFinder)

# Find ITK.

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

LINK_LIBRARIES(ITKBasicFilters ITKIO)

#CMake Optimizing
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(src)

SET(BoutonFinderSources src/bouton_params.cpp src/bouton_params.h src/bouton_finder.cpp src/bouton_finder.h src/bouton_finder_client.cpp src/amiraReader.cpp src/amiraReader.h src/typedefs.h src/basics.cpp src/basics.h src/typedefs_two.h src/utility.cpp src/utility.h src/brick_image_reader.cpp src/brick_image_reader.h)

#SET(NeuroMorphSources src/neuro_morph.cpp src/coordinates.cpp src/predecon.cpp src/streamer.cpp src/morphologyfilter.cpp src/objectlist.cpp src/compartmentlist.cpp src/compartment.cpp src/segmentation.cpp src/primitive.cxx src/skeletonization.cxx src/bloodvesselpattern.cpp src/image_reader.cpp src/hoc_streamer.cpp)


ADD_EXECUTABLE( BoutonFinder ${BoutonFinderSources} )

TARGET_LINK_LIBRARIES(BoutonFinder ITKCommon ITKIO)
SET_TARGET_PROPERTIES(BoutonFinder PROPERTIES LINK_FLAGS "-openmp")
