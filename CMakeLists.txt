find_package(OpenMesh REQUIRED)

add_executable(heatisoline HeatFermat3D.cpp Spirals.cpp SpiralDt.cpp)

target_include_directories(heatisoline PRIVATE ${OPENMESH_INCLUDE_DIR})

target_link_libraries(heatisoline PRIVATE
        igl::core
        igl::opengl
        igl::opengl_glfw
        ${OPENMESH_LIBRARIES} )