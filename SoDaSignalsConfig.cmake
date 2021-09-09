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

FIND_PACKAGE(FFTW3 REQUIRED)


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

