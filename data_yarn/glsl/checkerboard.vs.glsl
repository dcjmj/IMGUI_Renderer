#version 330 core

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;

in vec4 position;

vec3 constant_norm = vec3(0.0, 0.0, 1.0);

out VS_OUT
{
    vec2 tc;
	vec3 N;
	vec3 P;
} vs_out;

void main(void)
{
	const vec4 vertices[] =  vec4[](vec4( -2.0, -2.0, -0.1, 1.0),  
									vec4( -2.0,  2.0, -0.1, 1.0),   
									vec4(  2.0, -2.0, -0.1, 1.0),                                                        
                                    vec4(  2.0,  2.0, -0.1, 1.0));     

	gl_Position = proj_matrix * mv_matrix * vertices[gl_VertexID];
	vs_out.tc = vertices[gl_VertexID].xy;
	vs_out.N = normalize( constant_norm );  	
	vs_out.P = vertices[gl_VertexID].xyz;
}

