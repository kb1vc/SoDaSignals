include(FetchContent)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)

FetchContent_Declare(
  sodautils
  GIT_REPOSITORY https://github.com/kb1vc/SoDaUtils.git
  GIT_TAG v_3.0.0
)

FetchContent_MakeAvailable(sodautils)

