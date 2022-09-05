FIND_PACKAGE(SoDaUtils REQUIRED QUIET)


add_executable(FFTTest FFTTest.cxx)
target_link_libraries(FFTTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(FFTTiming FFTTiming.cxx)
target_link_libraries(FFTTiming sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FFTTiming PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FFTTiming PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(FilterTest FilterTest.cxx Checker.cxx)
target_link_libraries(FilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(FilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(FilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})

add_executable(OSFilterTest OSFilterTest.cxx Checker.cxx)
target_link_libraries(OSFilterTest sodasignals ${SoDaUtils_LIBRARIES})
target_link_directories(OSFilterTest PUBLIC ${SoDaUtils_LIBRARY_DIRS})
target_include_directories(OSFilterTest PUBLIC ${PROJECT_SOURCE_DIR}/include ${SoDaUtils_INCLUDE_DIRS})


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

add_test(NAME ReSamplerTest_1
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 48e3 --fsout 8e3)
set_tests_properties(ReSamplerTest_1 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_2
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 625e3 --fsout 48e3)
set_tests_properties(ReSamplerTest_2 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_3
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 120e3 --fsout 48e3)
set_tests_properties(ReSamplerTest_3 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_4
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 1.2e6 --fsout 48e3)
set_tests_properties(ReSamplerTest_4 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_5
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 1.25e6 --fsout 48e3)
set_tests_properties(ReSamplerTest_5 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_6
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsin 1.25e6 --fsout 48e3)
set_tests_properties(ReSamplerTest_6 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")


add_test(NAME ReSamplerTest_I1
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 48e3 --fsin 8e3)
set_tests_properties(ReSamplerTest_I1 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_I2
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 625e3 --fsin 48e3)
set_tests_properties(ReSamplerTest_I2 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_I3
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 120e3 --fsin 48e3)
set_tests_properties(ReSamplerTest_I3 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_I4
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 1.2e6 --fsin 48e3)
set_tests_properties(ReSamplerTest_I4 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_I5
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 1.25e6 --fsin 48e3)
set_tests_properties(ReSamplerTest_I5 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")

add_test(NAME ReSamplerTest_I6
  COMMAND $<TARGET_FILE:ReSamplerTest> --fsout 1.25e6 --fsin 48e3)
set_tests_properties(ReSamplerTest_I6 PROPERTIES
  PASS_REGULAR_EXPRESSION "PASSED")



