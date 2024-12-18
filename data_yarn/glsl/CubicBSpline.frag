#version 330 core

layout(location = 0, index = 0) out vec4 color;

void main()
{
	color.a = 1;
	color.rgb = vec3(0.68f,0.62f,0.18f);
}