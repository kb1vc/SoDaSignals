set(SIGNALS_SRCS
	FFT.cxx
	Periodogram.cxx
	Filter.cxx
	FilterSpec.cxx
	ReSampler.cxx
	NCO.cxx
	OSFilter.cxx
	Hilbert.cxx
	)


add_library(sodasignals STATIC ${SIGNALS_SRCS})
target_include_directories(sodasignals PUBLIC ${DEP_INCLS})
target_link_libraries(sodasignals ${DEP_LIBS})
target_link_directories(sodasignals PUBLIC ${FFTW3_LIBRARY_DIRS} ${FFTW3F_LIBRARY_DIRS})
set(SIGNALS_VERSION "${SoDaSignals_VERSION_MAJOR}.${SoDaSignals_VERSION_MINOR}.${SoDaSignals_VERSION_PATCH}")
set_target_properties(sodasignals PROPERTIES VERSION ${SIGNALS_VERSION})
install(TARGETS sodasignals DESTINATION ${CMAKE_INSTALL_LIBDIR})


