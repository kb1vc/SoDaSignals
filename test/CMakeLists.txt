FIND_PACKAGE(SoDaUtils REQUIRED QUIET)


add_executable(FFTTest FFTTest.cxx)
target_link_libraries(FFTTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(FFTTiming FFTTiming.cxx)
target_link_libraries(FFTTiming sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTiming PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTiming PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


add_executable(PeriodogramTest PeriodogramTest.cxx)
target_link_libraries(PeriodogramTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(PeriodogramTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(PeriodogramTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(ReSamplerTest ReSamplerTest.cxx)
target_link_libraries(ReSamplerTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(ReSamplerTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(ReSamplerTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(NCOTest NCOTest.cxx)
target_link_libraries(NCOTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(NCOTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(NCOTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


enable_testing()

add_test(NAME PeriodogramTest_1
  COMMAND $<TARGET_FILE:PeriodogramTest> --fsamp 48e3 -ftest 8e3 --psize 4096)
set_tests_properties(PeriodogramTest_1 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")
