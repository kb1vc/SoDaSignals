CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0) # this lets out really old distros like CentOS 6

IF(CMAKE_VERSION VERSION_GREATER 3.0.0)
  CMAKE_POLICY(SET CMP0048 NEW)
  SET(CMAKE_CXX_STANDARD 11)
ENDIF()

PROJECT(SoDaSignalExample)

# optimizations
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math")

# This is all you need to do if the SoDaFormat stuff
# was installed correctly. 
FIND_PACKAGE(SoDaSignals REQUIRED)

########### next target ###############


# Export header files that are required to build loadable modules

set(LIB_DIRS ${SoDaSignals_LIBRARY_DIRS})
set(LIBS ${SoDaSignals_LIBRARIES})
set(INCLUDE_DIRS ${SoDaSignals_INCLUDE_DIRS})

add_executable(FFTExample FFTExample.cxx)
target_link_directories(FFTExample PUBLIC ${LIB_DIRS})
target_link_libraries(FFTExample ${LIBS})
target_include_directories(FFTExample PUBLIC ${INCLUDE_DIRS})

add_executable(FilterExample FilterExample.cxx)
target_link_directories(FilterExample PUBLIC ${LIB_DIRS})
target_link_libraries(FilterExample ${LIBS})
target_include_directories(FilterExample PUBLIC ${INCLUDE_DIRS})


add_executable(FFTTiming FFTTiming.cxx)
target_link_directories(FFTTiming PUBLIC ${LIB_DIRS})
target_link_libraries(FFTTiming ${LIBS})
target_include_directories(FFTTiming PUBLIC ${INCLUDE_DIRS})

#add_executable(ReSamplerTest ReSamplerTest.cxx)
#target_link_directories(ReSamplerTest PUBLIC ${LIB_DIRS})
#target_link_libraries(ReSamplerTest ${LIBS})
#target_include_directories(ReSamplerTest PUBLIC ${INCLUDE_DIRS})



add_executable(FilterTest FilterTest.cxx)
target_link_directories(FilterTest PUBLIC ${LIB_DIRS})
target_link_libraries(FilterTest ${LIBS})
target_include_directories(FilterTest PUBLIC ${INCLUDE_DIRS})



