# Borrowed from https://wiki.o-ran-sc.org/display/ORAN/Packaging+Libraries+for+Linux+Distributions+with+CMake

IF( EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake" )
  include( InstallRequiredSystemLibraries )

  set( CPACK_SET_DESTDIR "off" )
  set( CPACK_PACKAGING_INSTALL_PREFIX "${install_root}" )
  set( CPACK_GENERATOR "DEB;RPM" )
 
  set( CPACK_PACKAGE_DESCRIPTION "SoDa package for basic signal processing functions (filters, resamplers, fft...)" )
  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "SoDa Signal Processing" )
  set( CPACK_PACKAGE_VENDOR "kb1vc" )
  set( CPACK_PACKAGE_CONTACT "https://github.com/kb1vc/" )
  set( CPACK_PACKAGE_VERSION_MAJOR "${SoDaSignals_VERSION_MAJOR}")
  set( CPACK_PACKAGE_VERSION_MINOR "${SoDaSignals_VERSION_MINOR}")
  set( CPACK_PACKAGE_VERSION_PATCH "${SoDaSignals_VERSION_PATCH}")
  set( CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}${spoiled_str}" )
  set( CPACK_SOURCE_PACKAGE_FILE_NAME
    "vric${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
 
  # omit if not required
  set( CPACK_DEBIAN_PACKAGE_DEPENDS "g++ >= 5.0.0, sodautils, fftw3-dev")
  set( CPACK_RPM_PACKAGE_REQUIRES "gcc-g++ >= 5.0.0, sodautils, fftw-devel")
 
  set( CPACK_DEBIAN_PACKAGE_PRIORITY "optional" )
  set( CPACK_DEBIAN_PACKAGE_SECTION "ric" )
  set( CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR} )
  set( CPACK_RPM_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR} )
 
  INCLUDE( CPack )
ENDIF()
