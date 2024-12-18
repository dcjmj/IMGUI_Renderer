#version 430 core

#define M_PI_2_INV	0.159154943091895335768883763372f

layout (lines) in; 
layout (triangle_strip, max_vertices = 256) out; 

in vec3 yarn_dir[];

uniform mat4 view_matrix;
uniform vec3 view_dir;
uniform float tube_width = 0.16f;
uniform float mDelta_0;
uniform float flyaway_texture_height;
uniform float width_scale;
uniform float texture_offset;
uniform mat4 shadow_matrix;
#define SHADOW_BIAS 0.0001f

in vec4 w_shadowPos[];
in vec4 w_pos[];
in float ply_shift[];
out vec4 pos;
out vec4 shadowPos;
out vec3 yarn_Dir;
out vec2 uv;
vec4 uv_[2];

vec3 GetWDif(vec3 sDir, vec3 view_dir_)
{
	vec3 bitangent = normalize(cross(normalize(sDir), view_dir_));
	//thickness here should include for hair fiber
	return bitangent * tube_width * width_scale;
}

void main(void)         
{
	float one_over_flyaway_texture_height = 1.0f / flyaway_texture_height;

	uv_[0]			= vec4( ply_shift[0] + texture_offset, ply_shift[0] + texture_offset, 1.0f - 12.0f*one_over_flyaway_texture_height, one_over_flyaway_texture_height );
	uv_[1]			= vec4( ply_shift[1] + texture_offset, ply_shift[1] + texture_offset, uv_[0].zw );
		
	vec3 view_dir0 = normalize(view_dir -  gl_in[0].gl_Position.xyz);
	vec3 view_dir1 = normalize(view_dir -  gl_in[1].gl_Position.xyz);

	vec3 width_0 = GetWDif( yarn_dir[0], view_dir0);
	vec3 width_1 = GetWDif( yarn_dir[1], view_dir1);

	///----------------------------------------------------------------------
	//	v1
	///----------------------------------------------------------------------
	gl_Position = gl_in[0].gl_Position;
	gl_Position.xyz += width_0 * 0.5;
	shadowPos = shadow_matrix * (gl_Position + vec4(view_dir0*tube_width*0.5, 0.0));
	shadowPos.z -= SHADOW_BIAS;
	gl_Position = view_matrix * gl_Position;
	gl_Position.z -= SHADOW_BIAS;
	yarn_Dir =  yarn_dir[0];
	uv = uv_[0].xz;

	pos = w_pos[0];

	EmitVertex();

	///----------------------------------------------------------------------
	//	v2
	///----------------------------------------------------------------------
	gl_Position = gl_in[0].gl_Position;
	gl_Position.xyz -= width_0 * 0.5;
	shadowPos = shadow_matrix * (gl_Position + vec4(view_dir0*tube_width*0.5, 0.0));
	shadowPos.z -= SHADOW_BIAS;
	gl_Position = view_matrix * gl_Position;
	gl_Position.z -= SHADOW_BIAS;
	uv	= uv_[0].yw;
	EmitVertex();

	///----------------------------------------------------------------------
	//	v3
	///----------------------------------------------------------------------
	gl_Position = gl_in[1].gl_Position;
	gl_Position.xyz += width_1 * 0.5;
	shadowPos = shadow_matrix * (gl_Position + vec4(view_dir1*tube_width*0.5, 0.0));
	shadowPos.z -= SHADOW_BIAS;
	gl_Position = view_matrix * gl_Position;
	gl_Position.z -= SHADOW_BIAS;
	yarn_Dir =  yarn_dir[1];
	uv	= uv_[1].xz;
	pos = w_pos[1];

	EmitVertex();

	///----------------------------------------------------------------------
	//	v4
	///----------------------------------------------------------------------
	gl_Position = gl_in[1].gl_Position;
	gl_Position.xyz -= width_1 * 0.5;
	shadowPos = shadow_matrix * (gl_Position + vec4(view_dir1*tube_width*0.5, 0.0));
	shadowPos.z -= SHADOW_BIAS;
	gl_Position = view_matrix * gl_Position;
	gl_Position.z -= SHADOW_BIAS;
	uv	= uv_[1].yw;
	EmitVertex();

	EndPrimitive();
}    