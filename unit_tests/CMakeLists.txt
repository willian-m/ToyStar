find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
ADD_DEFINITIONS(-lvec3)

if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)

#List of test executables
set(vec3_test_exe "vec3_test.exe")
set(particle_test_exe "particle_test.exe")

#Variables of libraries used in the test
#find_library(vec3_lib vec3 ${PROJECT_BINARY_DIR}/lib)
#message("vec3 is at: ${vec3_lib}")

#Configure build
add_executable("${vec3_test_exe}" vec3_test.cxx)
target_link_libraries (${vec3_test_exe} 
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       vec3
                       vec2
                       vec1
                       )

add_executable("${particle_test_exe}" particle_test.cxx)
target_link_libraries (${particle_test_exe}
                       particle
                       vec3
                       vec2
                       vec1
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${sph_math_test_exe}" sph_math_test.cxx)
target_link_libraries (${sph_math_test_exe}
                       sph_math
                       vec3
                       vec2
                       vec1
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${particle_system_test_exe}" particle_system_test.cxx)
target_link_libraries (${particle_system_test_exe}
                       sph_math
                       eos_polytropic
                       particle
                       particle_system 
                       cell
                       vec3
                       vec2
                       vec1
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${rk4_test_exe}" rk4_test.cxx)
target_link_libraries (${rk4_test_exe}
                       particle_system
                       particle
                       integrator_rk4
                       sph_math
                       eos_polytropic
                       vec3
                       vec2
                       vec1
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${eos_test_exe}" eos_test.cxx)
target_link_libraries (${eos_test_exe}
                       eos_polytropic
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )
#set_target_properties (${vec3_test_exe} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
#et_target_properties (${vec3_test_exe} PROPERTIES LINK_FLAGS "-lvec3")
                       

