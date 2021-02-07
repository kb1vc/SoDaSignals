set(SIGNALS_SRCS
	FFT.cxx
	ChebyshevWindow.cxx
	Filter.cxx
	)


FIND_PACKAGE(FFTW3 REQUIRED)
FIND_PACKAGE(FFTW3f REQUIRED)

set(DEP_LIBS ${FFTW3_LIBRARIES})
list(APPEND DEP_LIBS ${FFTW3f_LIBRARIES})

set(DEP_INCLS ${FFTW3_INCLUDE_DIRS})

message("DEP_INCLS = ${DEP_INCLS}")
message("DEP_LIBS = ${DEP_LIBS}")

add_library(sodasignals SHARED ${SIGNALS_SRCS})
target_include_directories(sodasignals PUBLIC ${DEP_INCLS})
target_link_libraries(sodasignals ${DEP_LIBS})
target_link_directories(sodasignals PUBLIC ${FFTW3_LIBRARY_DIRS} ${FFTW3f_LIBRARY_DIRS})
install(TARGETS sodasignals DESTINATION lib)




