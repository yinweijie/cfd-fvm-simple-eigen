function(cfd_resolve_eigen out_var)
  set(_candidate_paths)

  if(DEFINED EIGEN3_INCLUDE_DIR AND NOT EIGEN3_INCLUDE_DIR STREQUAL "")
    list(APPEND _candidate_paths "${EIGEN3_INCLUDE_DIR}")
  endif()

  list(APPEND _candidate_paths
    "${PROJECT_SOURCE_DIR}/third_party/eigen"
    "/usr/include/eigen3"
    "/usr/local/include/eigen3"
    "/opt/homebrew/include/eigen3"
    "/opt/local/include/eigen3"
  )

  foreach(_path IN LISTS _candidate_paths)
    if(EXISTS "${_path}/Eigen/Core")
      set(${out_var} "${_path}" PARENT_SCOPE)
      message(STATUS "Using Eigen headers from: ${_path}")
      return()
    endif()
  endforeach()

  string(JOIN "\n" _help_text
    "Eigen headers were not found."
    "Provide one of the following before configuring:"
    "  1. -DEIGEN3_INCLUDE_DIR=/path/to/eigen3"
    "  2. vendor Eigen under ${PROJECT_SOURCE_DIR}/third_party/eigen"
    "  3. install Eigen headers into a standard include path"
  )

  message(FATAL_ERROR "${_help_text}")
endfunction()
