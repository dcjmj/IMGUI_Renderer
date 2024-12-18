#version 430 core

layout (lines) in; 
layout (triangle_strip, max_vertices = 256) out; 

uniform float fiber_thickness = 0.008f;

out GEO_OUT
{
	vec3 fiber_dir;
	vec4 pos;
	vec3 yarn_offset;	// vector to yarn central line	(used for core)
	float is_core;		// 1.0 is core, 0.0 is regular
	vec2 uv;
	vec2 offset_2d;
	float fiber_thickness;
} geo_out;

uniform mat4 view_matrix;
uniform mat4 modelMatrix;
uniform float center_offset;

vec2 GetWDif( vec2 sPos1, vec2 sPos2 )
{
	vec2 sDir  = normalize(sPos2 - sPos1);
	vec2 sNorm = vec2( sDir.y, -sDir.x );
	return fiber_thickness*sNorm*1.5;
}

void main()
{
	vec2 width = GetWDif( (view_matrix * modelMatrix * (gl_in[0].gl_Position)).xy, (view_matrix * modelMatrix * (gl_in[1].gl_Position)).xy);
	width = vec2(0, 0.004);
	vec3 fiberDir = (view_matrix * modelMatrix * gl_in[0].gl_Position - view_matrix * modelMatrix * gl_in[1].gl_Position).xyz;

	fiberDir = normalize((gl_in[1].gl_Position - gl_in[0].gl_Position).xyz);
	///----------------------------------------------------------------------
	//	v1
	///----------------------------------------------------------------------

	gl_Position = view_matrix * modelMatrix * (gl_in[0].gl_Position);// + vec4(0.5 * width, 0, 0));
	gl_Position.xy += 0.5 * width;
	gl_Position.z += 0.2;
	geo_out.fiber_dir = fiberDir;
	geo_out.pos = gl_Position;
	geo_out.is_core = 0.0;
	geo_out.uv = vec2(0.05, 321.0f/326.0f);
	geo_out.yarn_offset = vec3(0.0, gl_in[0].gl_Position.y + center_offset, gl_in[0].gl_Position.z);
	geo_out.fiber_thickness = fiber_thickness;

	EmitVertex();

	///----------------------------------------------------------------------
	//	v2
	///----------------------------------------------------------------------
	
	gl_Position.xy -= width;
	geo_out.pos = gl_Position;
	geo_out.uv = vec2(0.95f, 321.0f/326.0f);
	
	EmitVertex();
	
	///----------------------------------------------------------------------
	//	v3
	///----------------------------------------------------------------------
	
	gl_Position = view_matrix * modelMatrix * (gl_in[1].gl_Position);// + vec4(0.5 * width, 0, 0));
	gl_Position.xy += 0.5 * width;
	gl_Position.z += 0.2;
	geo_out.fiber_dir = fiberDir;
	geo_out.pos = gl_Position;
	geo_out.uv = vec2(0.05, 321.0f/326.0f);
	geo_out.yarn_offset = vec3(0.0, gl_in[1].gl_Position.y + center_offset, gl_in[1].gl_Position.z);
	geo_out.fiber_thickness = fiber_thickness;

	EmitVertex();

	///----------------------------------------------------------------------
	//	v4
	///----------------------------------------------------------------------
	
	gl_Position.xy -= width;
	geo_out.pos = gl_Position;
	geo_out.uv = vec2(0.95f, 321.0f/326.0f);
	
	EmitVertex();

	EndPrimitive();
}
