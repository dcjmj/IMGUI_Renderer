#version 430 core

layout(location = 0, index = 0) out vec4 color;

in GEO_OUT
{
	vec3 fiber_dir;
	vec4 pos;
	vec3 yarn_offset;	// vector to yarn central line	(used for core)
	float is_core;		// 1.0 is core, 0.0 is regular
	vec2 uv;
	vec2 offset_2d;
	float fiber_thickness;
} fs_in;

vec3 N;

vec3 z = vec3(0,0,1);
float radius_ply = 0.035f;
float radius_fiber_min = 0.033f;

void main()
{	
	color.xyz = (fs_in.fiber_dir + vec3(1.0, 1.0, 1.0)) * 0.5f;
	color.w = 1.0;	// 0.05 is used to get the negative value, because GPU will not return a negative value
}