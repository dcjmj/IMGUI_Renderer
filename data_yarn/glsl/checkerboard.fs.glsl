#version 330 core

#define SPECULAR_GLOSS 1.0f

uniform sampler2D tex_object;
uniform vec3 light_pos;
in VS_OUT
{
    vec2 tc;
	vec3 N;
	vec3 P;
} fs_in;

out vec4 color;

void main(void)
{
	vec3 L = normalize( light_pos.xyz -  fs_in.P ); 
	float diffuseTerm = max(dot(fs_in.N, L), 0.0);

	/*
	vec3 V = -normalize(fs_in.pos.xyz/fs_in.pos.w);
	vec3 H = normalize(V + L);
	float s = 0;
	float sc = dot(N, H);
	if ( sc > 0 ) {
		s = pow( sc, SPECULAR_GLOSS );
	}
	*/

	vec4 textureColor = texture(tex_object, fs_in.tc);

    color =  textureColor * diffuseTerm + textureColor * 0.2;
	//color = vec4(gl_FragCoord.z);
}
