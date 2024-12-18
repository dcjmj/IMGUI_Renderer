#version 330 core

layout(location = 0) in vec3 inVertex;
layout(location = 2) in vec3 inNormal;

uniform mat4 view_matrix;
uniform mat4 modelMatrix;


out vec3 norm;
out vec3 pos;

void main()
{
	gl_Position = view_matrix * modelMatrix * vec4(inVertex,1);

	//gl_Position =  vec4(inVertex,1);
	norm = inNormal;
	pos = inVertex;
}
