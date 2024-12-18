#version 430 core

layout (lines) in; 
layout (triangle_strip, max_vertices = 256) out; 

uniform float fiber_thickness = 0.008f;

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
} geo_out;

uniform mat4 view_matrix;
uniform mat4 modelMatrix;
uniform float center_offset;
uniform float scale;

vec2 GetWDif( vec2 sPos1, vec2 sPos2 )
{
	vec2 sDir  = normalize(sPos2 - sPos1);
	vec2 sNorm = vec2( sDir.y, -sDir.x );
	return fiber_thickness*sNorm;
}

void main()
{
	vec2 width = GetWDif( (view_matrix * modelMatrix * (gl_in[0].gl_Position)).xy, (view_matrix * modelMatrix * (gl_in[1].gl_Position)).xy);
	
	geo_out.reverse = 1.0f;
	///----------------------------------------------------------------------
	//	v1
	///----------------------------------------------------------------------

	gl_Position = view_matrix * modelMatrix * (gl_in[0].gl_Position);// + vec4(0.5 * width, 0, 0));
	gl_Position.xy += 0.5 * width;
	geo_out.pos = gl_Position;
	geo_out.world_pos = gl_Position;
	geo_out.center_3d = gl_in[0].gl_Position.xyz;
	geo_out.uv = vec2(0.05, 321.0f/326.0f);
	geo_out.offset_3d = vec3(0.0, gl_in[0].gl_Position.y + center_offset, gl_in[0].gl_Position.z) * 1.5;
	geo_out.ply_dir_3d = vec3(-1,0,0);
	geo_out.fbr_dir_3d = vec3(-1,0,0);
	geo_out.fiber_thickness = fiber_thickness;

	EmitVertex();

	///----------------------------------------------------------------------
	//	v2
	///----------------------------------------------------------------------
	
	gl_Position.xy -= width;
	geo_out.pos = gl_Position;
	geo_out.world_pos = gl_Position;
	geo_out.uv = vec2(0.95f, 321.0f/326.0f);
	
	EmitVertex();
	
	///----------------------------------------------------------------------
	//	v3
	///----------------------------------------------------------------------
	
	gl_Position = view_matrix * modelMatrix * (gl_in[1].gl_Position);// + vec4(0.5 * width, 0, 0));
	gl_Position.xy += 0.5 * width;
	geo_out.pos = gl_Position;
	geo_out.world_pos = gl_Position;
	geo_out.center_3d = gl_in[1].gl_Position.xyz;
	geo_out.uv = vec2(0.05, 321.0f/326.0f);
	geo_out.offset_3d = vec3(0.0, gl_in[1].gl_Position.y + center_offset, gl_in[1].gl_Position.z) * 1.5;
	geo_out.ply_dir_3d = vec3(-1,0,0);
	geo_out.fbr_dir_3d = vec3(-1,0,0);
	geo_out.fiber_thickness = fiber_thickness;

	EmitVertex();

	///----------------------------------------------------------------------
	//	v4
	///----------------------------------------------------------------------
	
	gl_Position.xy -= width;
	geo_out.pos = gl_Position;
	geo_out.world_pos = gl_Position;
	geo_out.uv = vec2(0.95f, 321.0f/326.0f);
	
	EmitVertex();

	EndPrimitive();
}
