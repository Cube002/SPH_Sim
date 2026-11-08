cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

function(add_learnsph_app target source)
    add_executable(${target} ${source})
    add_dependencies(${target} learnSPH)
    target_link_libraries(${target} learnSPH)
endfunction()

# Assignment 1: neighborhood-search benchmark.
add_learnsph_app(assignment1_search_speed assignment1_search_speed.cpp)

# Assignment 1: simple gravity/neighborhood demo.
add_learnsph_app(assignment1_gravity_sim assignment1_gravity_sim.cpp)

# Assignment 2/3: boundary sampling and fluid-boundary density.
add_learnsph_app(assignment3_boundary_density assignment3_boundary_density.cpp)

# Assignment 4: pressure, friction and body accelerations.
add_learnsph_app(assignment4_forces_demo assignment4_forces_demo.cpp)

# Assignment 5: first no-boundary baseline test case.
add_learnsph_app(assignment5_test1_no_forces assignment5_test1_no_forces.cpp)

# Assignment 5a: incremental feature tests without boundaries.
add_learnsph_app(assignment5_a_i_pressure assignment5_a_i_pressure.cpp)
add_learnsph_app(assignment5_a_ii_viscosity assignment5_a_ii_viscosity.cpp)
add_learnsph_app(assignment5_a_iii_two_blocks assignment5_a_iii_two_blocks.cpp)
add_learnsph_app(assignment5_a_iv_gravity assignment5_a_iv_gravity.cpp)
add_learnsph_app(assignment5_a_v_floor assignment5_a_v_floor.cpp)
add_learnsph_app(assignment5_a_vi_boundary_viscosity assignment5_a_vi_boundary_viscosity.cpp)

# Assignment 5b: dam break.
add_learnsph_app(assignment5_b_dam_break assignment5_b_dam_break.cpp)
