cmake_minimum_required(VERSION 3.0.0)

project(gmp_test)

 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

set(CMAKE_CXX_STANDARD 11) # Tell CMake to compile with C++11


find_package(catkin REQUIRED COMPONENTS
  roscpp
  roslib
  gmp_lib
)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp
  DEPENDS
)

include_directories(
  include
  /usr/include/eigen3
  ${catkin_INCLUDE_DIRS}
)

set(SRC_LIST
  src/utils/utils.cpp
)

set(LIBS
  ${catkin_LIBRARIES}
  armadillo
  gmp_lib
)

# ================  Simulation  ===================

add_executable( nDoF_gmp_test src/simulation/nDoF_gmp_test.cpp ${SRC_LIST})
target_link_libraries(nDoF_gmp_test ${LIBS} )

add_executable( orient_gmp_test src/simulation/orient_gmp_test.cpp ${SRC_LIST})
target_link_libraries(orient_gmp_test ${LIBS} )

# ===============  Optimization  ===============

add_executable( pos_gmp_update src/optimization/update_gmp/pos_gmp_update.cpp ${SRC_LIST})
target_link_libraries(pos_gmp_update ${LIBS} )

add_executable( orient_gmp_update src/optimization/update_gmp/orient_gmp_update.cpp ${SRC_LIST})
target_link_libraries(orient_gmp_update ${LIBS} )

add_executable( constr_opt_test src/optimization/constr_opt/constr_opt_test.cpp ${SRC_LIST})
target_link_libraries(constr_opt_test ${LIBS} )

# ===================================================

#[[
add_executable(${PROJECT_NAME}
  include/${PROJECT_NAME}/main.cpp
)

add_dependencies(${PROJECT_NAME}
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  armadillo
  gmp_lib
)
]]

# =====================================
