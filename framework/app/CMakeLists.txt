cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    main.cpp)
    
add_executable(learnSPH_app ${SOURCE_FILES})
add_dependencies(learnSPH_app learnSPH)
target_link_libraries(learnSPH_app learnSPH)

# Gravity simulation
add_executable(gravitysim main_gravitysim.cpp)
add_dependencies(gravitysim learnSPH)
target_link_libraries(gravitysim learnSPH)

# Temporary target for neighborhood search speed tests
add_executable(searchSpeedTest searchSpeedTest.cpp)
add_dependencies(searchSpeedTest learnSPH)
target_link_libraries(searchSpeedTest learnSPH)
