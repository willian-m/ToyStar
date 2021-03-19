if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)

add_library(vec2 SHARED vec2.cxx)
set_target_properties (vec2 PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (vec2 PROPERTIES VERSION "0.1.0")
set_target_properties(vec2 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

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

add_library(eos_polytropic SHARED eos_polytropic.cxx)
set_target_properties (eos_polytropic PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (eos_polytropic PROPERTIES VERSION "0.1.0")
set_target_properties(eos_polytropic PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

add_library(particle_system SHARED particle_system.cxx)
set_target_properties (particle_system PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (particle_system PROPERTIES VERSION "0.1.0")
set_target_properties(particle_system PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

add_library(integrator_rk4 SHARED integrator_rk4.cxx)
set_target_properties (integrator_rk4 PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (integrator_rk4 PROPERTIES VERSION "0.1.0")
set_target_properties(integrator_rk4 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

add_library(cell SHARED cell.cxx)
set_target_properties (cell PROPERTIES COMPILE_FLAGS "${CompileFlags}")
set_target_properties (cell PROPERTIES VERSION "0.1.0")
set_target_properties(cell PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib/)

#install(TARGETS vec3 DESTINATION ${CMAKE_HOME_DIRECTORY})
