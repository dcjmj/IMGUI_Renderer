cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_language(CUDA)
project(render)
find_package(TBB CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(GLEW CONFIG REQUIRED)

file(GLOB_RECURSE SOURCES
  src/*.cpp
  imguiZMO/*.cpp
  imgui/backends/imgui_impl_glfw.cpp
  imgui/backends/imgui_impl_opengl3.cpp
  imgui/*.cpp
  )
file(GLOB_RECURSE HEADERS
  src/*.h
  imguiZMO/*.h
  imgui/*.h
  imgui/backends/*.h
  )
add_executable(ImGui_Render ${SOURCES})
target_include_directories(ImGui_Render PUBLIC src imgui imgui/backends imguiZMO)
target_link_libraries(ImGui_Render PUBLIC 
    TBB::tbb 
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    fmt::fmt
    GLEW::GLEW
    glfw
    glm::glm
    )
set_property(TARGET ImGui_Render PROPERTY CUDA_ARCHITECTURES 86)
set_target_properties(ImGui_Render PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_definitions(ImGui_Render PRIVATE ANTTWEAKBAR_DLL)
if(MSVC) 
    set_property(TARGET ImGui_Render PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data_yarn)
endif()

