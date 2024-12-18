#version 430 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 view_matrix;

flat out vec3 vertex_normal;
out vec3 posWorld;

void main()
{    
	gl_Position = view_matrix * vec4(position,1);
    vertex_normal = normal;
    posWorld = position;
}
