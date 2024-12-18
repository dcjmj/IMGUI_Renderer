#version 330 core

#define SHADOW_BIAS 0.003f

in vec4 pos;
in vec4 shadowPos;

uniform sampler2DShadow shadow_tex;
uniform mat4			shadow_matrix;

#if 0

layout(location = 0, index = 0) out vec4 color;

void main()
{
	fragmentdepth = gl_FragCoord.z;
}

#else

layout(location = 0, index = 0) out vec4 color;

void main()
{
	//vec4 shadowPos = shadow_matrix * pos;
	//float visibility = texture(shadow_tex, vec3(shadowPos.xy, (shadowPos.z)/shadowPos.w - SHADOW_BIAS));

   // color = vec4(gl_FragCoord.z);

   color.xyz = vec3(1.0, 0.0, 0.0);

   //color.xyz = vec3( visibility );

   color.xyz = vec3(textureProj(shadow_tex, shadowPos));

   color.w = 1.0;
}

#endif