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
add_executable("${vec3_test_exe}" vec3_test)
target_link_libraries (${vec3_test_exe} 
                       vec3
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${particle_test_exe}" particle_test)
target_link_libraries (${particle_test_exe} 
                       vec3
                       particle
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${sph_math_test_exe}" sph_math_test)
target_link_libraries (${sph_math_test_exe} 
                       vec3
                       sph_math
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

add_executable("${particle_system_test_exe}" particle_system_test)
target_link_libraries (${particle_system_test_exe}
                       sph_math
                       eos_polytropic
                       particle
                       particle_system 
                       vec3
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       )

#set_target_properties (${vec3_test_exe} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
#et_target_properties (${vec3_test_exe} PROPERTIES LINK_FLAGS "-lvec3")
                       

