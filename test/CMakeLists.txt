FIND_PACKAGE(SoDaUtils  QUIET)


add_executable(CircularBufferTest CircularBufferTest.cxx)
if(SoDaUtils_FOUND) 
add_executable(ReSamplerTest ReSamplerTest.cxx)
target_link_libraries(ReSamplerTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(ReSamplerTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(ReSamplerTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(ReSamplerNoise ReSamplerNoise.cxx)
target_link_libraries(ReSamplerNoise sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(ReSamplerNoise PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(ReSamplerNoise PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


add_executable(FilterTest FilterTest.cxx)
target_link_libraries(FilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(FFTFilterTest FFTFilterTest.cxx)
target_link_libraries(FFTFilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTFilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTFilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


add_executable(FFTTest FFTTest.cxx)
target_link_libraries(FFTTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})




add_executable(ReSamplerSweep ReSamplerSweep.cxx)
target_link_libraries(ReSamplerSweep sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(ReSamplerSweep PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(ReSamplerSweep PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})
ENDIF()

add_executable(FFT_Timing FFT_Timing.cxx)
target_link_libraries(FFT_Timing sodasignals)
target_include_directories(FFT_Timing PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(FFT_Timing2 FFT_Timing2.cxx)
target_link_libraries(FFT_Timing2 sodasignals)
target_include_directories(FFT_Timing2 PUBLIC ${PROJECT_SOURCE_DIR}/include)
