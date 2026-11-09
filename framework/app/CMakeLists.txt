cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

function(add_learnsph_app target source)
    add_executable(${target} ${source})
    add_dependencies(${target} learnSPH)
    target_link_libraries(${target} learnSPH)
endfunction()

set(SHEET1_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Sheet1)
set(SHEET2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Sheet2)

# Sheet 1: neighborhood-search benchmark.
add_learnsph_app(assignment1_search_speed ${SHEET1_DIR}/assignment1_search_speed.cpp)

# Sheet 1: simple gravity/neighborhood demo.
add_learnsph_app(assignment1_gravity_sim ${SHEET1_DIR}/assignment1_gravity_sim.cpp)

# Sheet 2, Assignment 1: fluid-box density test case.
add_learnsph_app(assignment1_fluid_density ${SHEET2_DIR}/assignment1_fluid_density.cpp)

# Sheet 2, Assignment 2: boundary sampling and volume/mass correction.
add_learnsph_app(assignment2_boundary_sampling ${SHEET2_DIR}/assignment2_boundary_sampling.cpp)

# Sheet 2, Assignment 3: fluid-box density including a boundary floor.
add_learnsph_app(assignment3_fluid_boundary_density ${SHEET2_DIR}/assignment3_fluid_boundary_density.cpp)

# Sheet 2, Assignment 2/3: boundary sampling and fluid-boundary density.
add_learnsph_app(assignment3_boundary_density ${SHEET2_DIR}/assignment3_boundary_density.cpp)

# Sheet 2, Assignment 4: pressure, friction and body accelerations.
add_learnsph_app(assignment4_forces_demo ${SHEET2_DIR}/assignment4_forces_demo.cpp)

# Sheet 2, Assignment 5: first no-boundary baseline test case.
add_learnsph_app(assignment5_test1_no_forces ${SHEET2_DIR}/assignment5_test1_no_forces.cpp)

# Sheet 2, Assignment 5a: incremental feature tests without boundaries.
add_learnsph_app(assignment5_a_i_pressure ${SHEET2_DIR}/assignment5_a_i_pressure.cpp)
add_learnsph_app(assignment5_a_ii_viscosity ${SHEET2_DIR}/assignment5_a_ii_viscosity.cpp)
add_learnsph_app(assignment5_a_iii_two_blocks ${SHEET2_DIR}/assignment5_a_iii_two_blocks.cpp)
add_learnsph_app(assignment5_a_iv_gravity ${SHEET2_DIR}/assignment5_a_iv_gravity.cpp)
add_learnsph_app(assignment5_a_v_floor ${SHEET2_DIR}/assignment5_a_v_floor.cpp)
add_learnsph_app(assignment5_a_vi_boundary_viscosity ${SHEET2_DIR}/assignment5_a_vi_boundary_viscosity.cpp)

# Sheet 2, Assignment 5b: dam break.
add_learnsph_app(assignment5_b_dam_break ${SHEET2_DIR}/assignment5_b_dam_break.cpp)

# Sheet 2, Assignment 5c: large showcase scene.
add_learnsph_app(assignment5_c_obstacle_dam_break ${SHEET2_DIR}/assignment5_c_obstacle_dam_break.cpp)
