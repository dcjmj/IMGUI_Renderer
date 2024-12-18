#version 330 core

layout(location = 0) in vec2 inVertex;

void main()
{
	gl_Position = vec4(inVertex,0,0);
}
