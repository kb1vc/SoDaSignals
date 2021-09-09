#
#  Copyright (c) 2019 Matthew H. Reilly (kb1vc)
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#
#  Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#  Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in
#  the documentation and/or other materials provided with the
#  distribution.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

set(SoDaSignals_FOUND TRUE)

INCLUDE(FindPackageHandleStandardArgs)
find_package(PkgConfig)

set(SoDaSignals_INCLUDE_HINTS)
set(SoDaSignals_LIB_HINTS)



find_path(SoDaSignals_pre_INCLUDE_DIRS
  NAMES SoDa/Filter.hxx
  HINTS ${SoDaSignals_INCLUDE_HINTS}
  PATHS /usr/local/include
        /usr/include
	/opt/local/include
)

find_library(SoDaSignals_pre_LIBRARIES
  NAMES sodasignals
  HINTS ${SoDaSignals_LIB_HINTS}
  PATHS /usr/local/lib
       /usr/lib
       /usr/local/lib64
       /usr/lib64
       /opt/local/lib
)

# We need to use our own FFTW finder because the FFTW installation is
# broken on some Fedora releases... (It is missing a CMake Config file.)

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

SET(SoDaSignals_INCLUDE_DIRS ${SoDaSignals_pre_INCLUDE_DIRS})
SET(SoDaSignals_LIBRARIES ${SoDaSignals_pre_LIBRARIES})

LIST(APPEND SoDaSignals_INCLUDE_DIRS ${FFTW3_INCLUDE_DIRS})
LIST(APPEND SoDaSignals_LIBRARIES ${FFTW3_LIBRARIES})

if(SoDaSignals_pre_INCLUDE_DIRS AND SoDaSignals_pre_LIBRARIES)
  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(SoDaSignals DEFAULT_MSG SoDaSignals_LIBRARIES SoDaSignals_INCLUDE_DIRS)
  mark_as_advanced(SoDaSignals_LIBRARIES SoDaSignals_INCLUDE_DIRS)
else()
  message(FATAL_ERROR "SoDaSignals lib is required, but not found.")
endif()

