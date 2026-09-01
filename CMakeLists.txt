cmake_minimum_required(VERSION 3.16)
project(MeshSimProject LANGUAGES CXX)

# Options
option(BUILD_MESHSIM_CLI "Build the meshsim_cli executable" ON)
option(BUILD_BUFFER_CERTIFY "Build the buffer_certify tool (separate from MeshSim)" ON)
option(BUILD_RECHECK "Build the recheck_candidates tool" ON)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Reasonable defaults
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# -------- MeshSim library --------
add_library(meshsim STATIC
  MeshSim.cpp
)
target_include_directories(meshsim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(meshsim PUBLIC cxx_std_17)

# -------- meshsim_cli executable --------
if(BUILD_MESHSIM_CLI)
  add_executable(meshsim_cli
    MeshSim_cli.cpp
  )
  target_link_libraries(meshsim_cli PRIVATE meshsim)
endif()

# -------- buffer_certify tool (separate) --------
# Rename the source file below to match your actual filename in the repo:
#   buffer_hypothesis_certify_3phases_iut.cpp
#
if(BUILD_BUFFER_CERTIFY)
  add_executable(buffer_certify
    buffer_hypothesis_certify_3phases_iut.cpp
  )
  target_link_libraries(buffer_certify PRIVATE meshsim)
endif()

# -------- recheck_candidates tool --------
# Re-simulates the leading Phase-3 competitors at a long common horizon.
if(BUILD_RECHECK)
  find_package(Threads REQUIRED)
  add_executable(recheck_candidates
    recheck_candidates.cpp
  )
  target_link_libraries(recheck_candidates PRIVATE meshsim Threads::Threads)
endif()

# -------- Compiler warnings (optional, safe) --------
if(MSVC)
  target_compile_options(meshsim PRIVATE /W3)
  if(BUILD_MESHSIM_CLI)
    target_compile_options(meshsim_cli PRIVATE /W3)
  endif()
  if(BUILD_BUFFER_CERTIFY)
    target_compile_options(buffer_certify PRIVATE /W3)
  endif()
  if(BUILD_RECHECK)
    target_compile_options(recheck_candidates PRIVATE /W3)
  endif()
else()
  target_compile_options(meshsim PRIVATE -Wall -Wextra -Wpedantic)
  if(BUILD_MESHSIM_CLI)
    target_compile_options(meshsim_cli PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  if(BUILD_BUFFER_CERTIFY)
    target_compile_options(buffer_certify PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  if(BUILD_RECHECK)
    target_compile_options(recheck_candidates PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()