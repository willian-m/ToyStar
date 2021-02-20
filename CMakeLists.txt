cmake_minimum_required(VERSION 3.0.0)
project(ToyStar VERSION 0.1.0 LANGUAGES CXX)

#configure a header file to pass the version number to the source code:
#configure_file(TutorialConfig.h.in TutorialConfig.h)

# default build type: Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_BUILD_TYPE "Debug")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_HOME_DIRECTORY}/"
        CACHE PATH "Install path prefix, prepended onto install directories."
        FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

#Include folders
include_directories(${PROJECT_SOURCE_DIR}/include)
#Library folder
link_directories(${PROJECT_BINARY_DIR}/lib)
#Executable folder
set(vec3_test_exe "vec3_test.exe")
set(particle_test_exe "particle_test.exe")
set(sph_math_test_exe "sph_math_test.exe")
set(particle_system_test_exe "particle_system_test.exe")
set(rk4_test_exe "rk4_test.exe")

string(APPEND CMAKE_CXX_FLAGS " -Wall")

add_subdirectory (src)
add_subdirectory (unit_tests)

enable_testing()
add_test(Vec3 ${PROJECT_BINARY_DIR}/unit_tests/${vec3_test_exe})
add_test(Particle ${PROJECT_BINARY_DIR}/unit_tests/${particle_test_exe})
add_test(SPHMath ${PROJECT_BINARY_DIR}/unit_tests/${sph_math_test_exe})
add_test(ParticleSystem ${PROJECT_BINARY_DIR}/unit_tests/${particle_system_test_exe})
add_test(RungeKutta4Integrator ${PROJECT_BINARY_DIR}/unit_tests/${rk4_test_exe})