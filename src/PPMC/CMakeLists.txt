# PPMC CMake script
# Adrien Maglo

cmake_minimum_required(VERSION 2.6)

project(PPMC)

# Find CGAL, OpenGL and GLUT
find_package(CGAL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include(${CGAL_USE_FILE})

add_executable(ppmc rangeCoder/rangecod.c rangeCoder/qsmodel.c main.cpp mymesh.cpp mymeshComp.cpp mymeshCompTests.cpp mymeshDecomp.cpp mymeshAdaptiveQuantization.cpp mymeshLifting.cpp mymeshUtils.cpp frenetRotation.cpp mymeshIO.cpp)

set_target_properties(ppmc PROPERTIES COMPILE_FLAGS -frounding-math)

target_link_libraries(ppmc ${CGAL_LIBRARIES} ${OPENGL_LIBRARY} ${GLUT_LIBRARY} m)

