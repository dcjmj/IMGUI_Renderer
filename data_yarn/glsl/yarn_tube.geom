#version 430 core

#define M_PI_2_INV	0.159154943091895335768883763372f

layout (lines) in; 
layout (triangle_strip, max_vertices = 256) out; 

in vec3 yarn_dir[];

uniform mat4 view_matrix;
uniform vec3 light_dir;
uniform float tube_width = 0.16f;

in vec4 w_shadowPos[];
in vec4 w_pos[];
out vec4 pos;
out vec4 shadowPos;

vec3 GetWDif(vec3 sDir)
{
	vec3 bitangent = normalize(cross(normalize(sDir), light_dir));
	return bitangent * tube_width;
}

void main(void)         
{
	vec3 width_0 = GetWDif( yarn_dir[0] );
	vec3 width_1 = GetWDif( yarn_dir[1] );

	///----------------------------------------------------------------------
	//	v1
	///----------------------------------------------------------------------
	gl_Position = gl_in[0].gl_Position;
	gl_Position.xyz += width_0 * 0.5;
	gl_Position = view_matrix * gl_Position;

	pos = w_pos[0];
	shadowPos = w_shadowPos[0];

	EmitVertex();

	///----------------------------------------------------------------------
	//	v2
	///----------------------------------------------------------------------
	gl_Position = gl_in[0].gl_Position;
	gl_Position.xyz -= width_0 * 0.5;
	gl_Position = view_matrix * gl_Position;
	EmitVertex();

	///----------------------------------------------------------------------
	//	v3
	///----------------------------------------------------------------------
	gl_Position = gl_in[1].gl_Position;
	gl_Position.xyz += width_1 * 0.5;
	gl_Position = view_matrix * gl_Position;

	pos = w_pos[1];
	shadowPos = w_shadowPos[1];

	EmitVertex();

	///----------------------------------------------------------------------
	//	v4
	///----------------------------------------------------------------------
	gl_Position = gl_in[1].gl_Position;
	gl_Position.xyz -= width_1 * 0.5;
	gl_Position = view_matrix * gl_Position;
	EmitVertex();

	EndPrimitive();
}    