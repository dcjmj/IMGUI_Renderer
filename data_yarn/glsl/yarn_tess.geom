#version 430 core

//#define HIGH_PERFORMANCE

layout (lines) in; 
layout (triangle_strip, max_vertices = 256) out; 

in TES_OUT
{	
	bool is_core;
	bool fly_away;
	vec3 core_dir;
	float ply_shift;
	float fiber_thickness;
	vec3 offset_3d;
	vec3 ply_dir_3d;
	vec3 fbr_dir_3d;
	vec3 center_3d;				// world space center
	vec4 world_pos;
	vec3 lod_color;
	vec3 up;
	float curvature;
} geo_in[];

out GEO_OUT
{
	vec4 pos;
	vec2 uv;
	float fiber_thickness;
	vec3 offset_3d;
	vec3 ply_dir_3d;
	vec3 fbr_dir_3d;
	vec3 center_3d;	
	float reverse;
	vec4 world_pos;
	vec3 lod_color;
	vec3 up;
	float curvature;
} geo_out;
uniform float mDelta_1;

uniform float fiber_thickness;
uniform float ellipse_long;
uniform float ellipse_short;
//uniform float scale;
uniform float core_texture_height;

uniform bool use_regular_fiber;

uniform mat4 view_matrix;
uniform vec3 view_pos;
uniform float mDelta_0;

vec3 GetWDif(vec3 sDir, vec3 view_dir)
{
	vec3 bitangent = normalize(cross(normalize(sDir), view_dir));

	return bitangent;
}

float coreCentralHighlightThreshold = 0.5f;
vec3 width[2];
vec4 uv[2];
vec4 pos[2];

void main(void)         
{
	bool reverse_z = view_matrix[0][0] > 0.0;	//(dot( mat3(view_matrix) * vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f) ) > 0.0);

#ifndef HIGH_PERFORMANCE	
	if(use_regular_fiber || geo_in[1].is_core)	// DISABLE HERE TO BE FASTER !!!
#endif
	{
		if (!geo_in[0].fly_away && !geo_in[1].fly_away)	
		{
			geo_out.reverse		= reverse_z ? 1.0 : -1.0;
			geo_out.lod_color	= geo_in[0].lod_color;
			

			pos[0]				= gl_in[0].gl_Position;
			pos[1]				= gl_in[1].gl_Position;

			vec3 view_dir0 = normalize(view_pos - pos[0].xyz);
			vec3 view_dir1 = normalize(view_pos - pos[1].xyz);

			width[0]		= GetWDif( geo_in[0].core_dir, view_dir0 ) * geo_in[0].fiber_thickness;
			width[1]		= GetWDif( geo_in[1].core_dir, view_dir1 ) * geo_in[1].fiber_thickness;

			float one_over_core_texture_height = 1.0f / core_texture_height;

			if ( geo_in[0].is_core ) {
				uv[0]			= vec4( geo_in[0].ply_shift + mDelta_0, geo_in[0].ply_shift + mDelta_0, 1.0f - 11.0f*one_over_core_texture_height, one_over_core_texture_height );
				uv[1]			= vec4( geo_in[1].ply_shift + mDelta_0, geo_in[1].ply_shift + mDelta_0, uv[0].zw );
			} else {
				float f = 1.0f - 5.5f * one_over_core_texture_height;
				uv[0]			= vec4( reverse_z ? vec2(0.95,0):vec2(0,0.95), f, f );
				uv[1]			= uv[0];
			}

			///----------------------------------------------------------------------
			//	v1
			///----------------------------------------------------------------------
			gl_Position			= pos[0];
			gl_Position.xyz		+= 0.5 * width[0];
			geo_out.world_pos	=  vec4(gl_Position.xyz, 1.0);
			gl_Position = view_matrix * gl_Position;
	
			geo_out.pos			= pos[0];
			geo_out.uv			= uv[0].xz;
			geo_out.offset_3d	= geo_in[0].offset_3d;
			geo_out.ply_dir_3d	= geo_in[0].ply_dir_3d;
			geo_out.fbr_dir_3d	= geo_in[0].fbr_dir_3d;
			geo_out.center_3d	= geo_in[0].center_3d;
			geo_out.fiber_thickness = geo_in[0].fiber_thickness;
			geo_out.up			= geo_in[0].up;
			geo_out.curvature			= geo_in[0].curvature;
			EmitVertex();

			///----------------------------------------------------------------------
			//	v2
			///----------------------------------------------------------------------
			gl_Position			= pos[0];
			gl_Position.xyz		-= 0.5 * width[0];
			geo_out.world_pos	=  vec4(gl_Position.xyz, 1.0);
			gl_Position = view_matrix * gl_Position;
			geo_out.uv			= uv[0].yw;
			EmitVertex();

			///----------------------------------------------------------------------
			//	v3
			///----------------------------------------------------------------------
			gl_Position			= pos[1];
			gl_Position.xyz		+= 0.5 * width[1];
			geo_out.world_pos	=  vec4(gl_Position.xyz, 1.0);
			gl_Position = view_matrix * gl_Position;

			geo_out.pos			= pos[1];
			geo_out.uv			= uv[1].xz;
			geo_out.offset_3d	= geo_in[1].offset_3d;
			geo_out.ply_dir_3d	= geo_in[1].ply_dir_3d;
			geo_out.fbr_dir_3d	= geo_in[1].fbr_dir_3d;
			geo_out.center_3d	= geo_in[1].center_3d;
			geo_out.fiber_thickness = geo_in[1].fiber_thickness;
			geo_out.up			= geo_in[1].up;
			geo_out.curvature			= geo_in[1].curvature;
			EmitVertex();

			///----------------------------------------------------------------------
			//	v4
			///----------------------------------------------------------------------
			gl_Position			= pos[1];
			gl_Position.xyz		-= 0.5 * width[1];
			geo_out.world_pos	=  vec4(gl_Position.xyz, 1.0);
			gl_Position = view_matrix * gl_Position;
			geo_out.uv			= uv[1].yw;
			EmitVertex();

			EndPrimitive();
		}	
	}
}    