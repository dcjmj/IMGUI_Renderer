#version 330 core

#define HAIR_SHADOW_LAYER_SIZE 0.05f
#define HAIR_SHADOW_OPACITY 0.03f

// Ouput data
layout(location = 0, index = 0) out vec4 color;

in vec4 pos;

uniform sampler2D hairShadow;

void main()
{
	vec2 tc = pos.xy*0.5f + vec2(0.5f,0.5f);
	float hairZ = texture(hairShadow, tc).r;
	float z = pos.z*0.5f + 0.5f;
	float zDif = abs(hairZ-z);
	int layer = min( 3, int( zDif * (1.0f/HAIR_SHADOW_LAYER_SIZE) ) );

	vec4 shad[4] = vec4[]( vec4(1,1,1,1), vec4(0,1,1,1), vec4(0,0,1,1), vec4(0,0,0,1) );

	color = shad[ layer ] * HAIR_SHADOW_OPACITY;
}