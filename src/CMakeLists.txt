if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)

add_library(vec3 SHARED vec3.cxx)
set_target_properties (vec3 PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (vec3 PROPERTIES VERSION "0.1.0")
set_target_properties(vec3 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

add_library(particle SHARED particle.cxx)
set_target_properties (particle PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (particle PROPERTIES VERSION "0.1.0")
set_target_properties(particle PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

add_library(sph_math SHARED sph_math.cxx)
set_target_properties (sph_math PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (sph_math PROPERTIES VERSION "0.1.0")
set_target_properties(sph_math PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)


#install(TARGETS vec3 DESTINATION ${CMAKE_HOME_DIRECTORY})
