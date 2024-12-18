#version 430 core

layout(location = 0, index = 0) out vec4 color;


in float ran;
in vec3 yarnDir;
in vec4 pos;
in vec3 norm;
in vec4 shadow_coord;
in float distance2centralLine;
in vec2 texture_coord;


void main()
{
	color.a = 1;
	color.rgb = vec3(0.68f,0.62f,0.18f);
}