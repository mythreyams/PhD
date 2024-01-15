# This project is designed to be built outside the Insight source tree.
PROJECT(ErrorCalculator)


#CMake Optimizing
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_C_FLAGS_RELEASE "-O2 -openmp")
SET(CMAKE_VERBOSE_MAKEFILE "TRUE")

#Sources
INCLUDE_DIRECTORIES(src)

SET(ErrorCalculatorSources  src/bouton_finder_client.cpp src/result_analyzer.cpp src/result_analyzer.h)

#SET(NeuroMorphSources src/neuro_morph.cpp src/coordinates.cpp src/predecon.cpp src/streamer.cpp src/morphologyfilter.cpp src/objectlist.cpp src/compartmentlist.cpp src/compartment.cpp src/segmentation.cpp src/primitive.cxx src/skeletonization.cxx src/bloodvesselpattern.cpp src/image_reader.cpp src/hoc_streamer.cpp)


ADD_EXECUTABLE( ErrorCalculator ${ErrorCalculatorSources} )


SET_TARGET_PROPERTIES(ErrorCalculator PROPERTIES LINK_FLAGS "-openmp")
