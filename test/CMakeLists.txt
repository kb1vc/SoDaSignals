FIND_PACKAGE(SoDaUtils REQUIRED)

add_executable(CircularBufferTest CircularBufferTest.cxx)

#add_executable(ReSamplerTest ReSamplerTest.cxx)
#target_link_libraries(ReSamplerTest sodasignals ${SoDaUtils_LIBRARIES})
#target_link_directories(ReSamplerTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
#target_include_directories(ReSamplerTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


add_executable(FilterTest FilterTest.cxx)
target_link_libraries(FilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(FIRFilterTest FIRFilterTest.cxx)
target_link_libraries(FIRFilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FIRFilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FIRFilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})



add_executable(FFTTest FFTTest.cxx)
target_link_libraries(FFTTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})




