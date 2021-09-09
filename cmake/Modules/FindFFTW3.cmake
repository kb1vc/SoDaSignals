# The Fedora (34?) FFTW package is missing some CMake support.
# In particular the FFTW3LibraryDepends.cmake file is missing.
# There are various hacks to fix this, but I think we'll just bypass
# the whole thing. 

FIND_PATH(
    FFTW3_INCLUDE_DIRS
    NAMES fftw3.h 
    HINTS $ENV{FFTW3_DIR}/include
    PATHS /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    FFTW3F_libs
    NAMES fftw3f libfftw3f
    HINTS $ENV{FFTW3_DIR}/lib
    PATHS /usr/local/lib
          /usr/lib
          /usr/lib64
)

FIND_LIBRARY(
    FFTW3_libs
    NAMES fftw3 libfftw3
    HINTS $ENV{FFTW3_DIR}/lib
    PATHS /usr/local/lib
          /usr/lib
          /usr/lib64
)

set(FFTW3_LIBRARIES ${FFTW3_libs} ${FFTW3F_libs})
	
INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIRS)

MARK_AS_ADVANCED(FFTW3_LIBRARIES FFTW3_INCLUDE_DIRS)

