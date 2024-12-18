#version 330 core

#define LIGHT_INTENSITY 1.0f
#define AMBIENT_LIGHT vec3(0.1f,0.12f,0.14f)
#define DIFFUSE_COLOR vec3(0.68f,0.62f,0.18f)
#define SPECULAR_COLOR (vec3(1,1,1)*0.7f)
#define SPECULAR_GLOSS 40.0f

layout(location = 0, index = 0) out vec4 color;

in vec3 normal;
in vec4 pos;

uniform vec3 lightDir;

void main()
{
	color.a = 1;

	vec3 dLight = AMBIENT_LIGHT;
	float sLight = 0;
	vec3 N = normalize(normal);
	vec3 L = lightDir;
	float c = dot(N, L);
	if ( c > 0 ) {
		dLight = LIGHT_INTENSITY * vec3(c,c,c);
		vec3 V = normalize(pos.xyz/pos.w);
		vec3 H = normalize(V + L);
		float sc = dot(N, H);
		if ( sc > 0 ) {
			sLight = LIGHT_INTENSITY * c * pow( sc, SPECULAR_GLOSS );
		}
	}
	color.rgb = DIFFUSE_COLOR * dLight + sLight*SPECULAR_COLOR;
}