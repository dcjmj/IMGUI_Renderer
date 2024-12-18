#version 430 core

layout(location = 0, index = 0) out vec4 color;

in GEO_OUT
{
	vec4 pos;
	vec2 uv;
	vec2 offset_2d;
	float fiber_thickness;
} fs_in;

uniform vec3 lightPos;
uniform sampler2D core_tex;
uniform bool use_ao;
uniform bool use_diffuse;
uniform bool use_regular_fiber;
uniform float scale;

vec3 V, N, L;
vec3 light_dir = vec3(-1.0, -1.0, 1.0);
vec3 z = vec3(0,0,1);
vec3 defaultColor = vec3(0.2f, 0.8f,0.1f);
float max_fiber_radius = 0.07f;

void main()
{
	// light direction
    L = normalize(light_dir);

	// fetch texture
	vec4 tex_value = texture(core_tex, fs_in.uv);

	N = tex_value.xyz * 2.0f - vec3(1.0, 1.0, 1.0);

	if (tex_value.w < 0.5)	discard;

	// diffuse
	float diffuse = max(dot(N, L), 0.0);

	///----------------------------------------------------------------------		
	//	ambient occlusion
	//	texture_norm.w			range [-1, 1]
	//	dist2FiberCenter		range [-1, 1]
	///----------------------------------------------------------------------
		
	float dist2FiberCenter = (0.5f - fs_in.uv.y) * 2.0f;	
	vec2 height = vec2(dist2FiberCenter, tex_value.w * 4.0f - 3.0f) * fs_in.fiber_thickness * 0.5f;

	float ao = length( (fs_in.offset_2d + height) / 2.0) / max_fiber_radius * 1.3;

	ao = min(max(ao * ao * 5.0f , 0.2f), 1.0f);

	if(!use_ao)	ao = 1.0f;
	if(!use_diffuse) diffuse = 1.0f;

	// finalize
	color.xyz = diffuse * defaultColor * ao;
	color.w = 1.0;
}