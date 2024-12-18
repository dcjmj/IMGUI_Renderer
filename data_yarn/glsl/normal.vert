#version 330 core

#define LIGHT_POS vec3(100,200,1000)

layout(location = 0) in vec3 inVertex;
layout(location = 2) in vec3 inNormal;

out vec3 normal;
out vec4 pos;

uniform mat4 cameraMatrix;
uniform mat3 cameraNormalMatrix;

void main()
{
	gl_Position = cameraMatrix * vec4(inVertex,1);
	normal      = cameraNormalMatrix * inNormal;
	pos         = gl_Position;
}
