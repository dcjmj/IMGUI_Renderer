-- xmake.lua
add_rules("mode.debug")
set_languages("c++17")
set_runtimes("MT")

add_requires("cuda", {system=true})
add_requires("eigen >=3.4.0")
add_requires("fmt =8.1.1")
add_requires("boost =1.78.0")
add_requires("nlohmann_json =3.10.5")
add_requires("tbb")
add_requires("opengl")
add_requires("glfw")
add_requires("glew")
add_requires("glu")
add_requires("glm")

target("ImGui_Render")
    set_kind("binary")
    add_packages("cuda")
    add_packages("eigen", "fmt", "boost", "nlohmann_json", "tbb")
    add_packages("opengl", "glfw", "glew", "glu", "glm")

    add_includedirs("src")
    add_includedirs("imgui", "imgui/backends", "imguiZMO")
    add_files("src/*.cpp")
    add_files("imguiZMO/*.cpp")
    add_files("imgui/backends/imgui_impl_glfw.cpp", "imgui/backends/imgui_impl_opengl3.cpp")
    add_files("imgui/*.cpp")
    
    add_defines("ANTTWEAKBAR_DLL")
    if is_plat("windows") then
        add_cuflags("-rdc=true", {force = true})
        --add_ldflags("/utf-8", {force = true})
        --add_cxflags("/utf8")
    end