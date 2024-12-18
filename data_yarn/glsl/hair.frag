#version 330 core

#define LIGHT_INTENSITY 1.0f
#define AMBIENT_LIGHT (vec3(0.12f,0.12f,0.14f))
#define DIFFUSE_COLOR ( vec3(93,40,12)/255.0f*0.5f )	// brown
#define SPECULAR_COLOR (vec3(1,1,1)*0.6f)
#define SPECULAR_GLOSS 80.0f

#define SHADOW_BIAS 0.005f

#define HAIR_SHADOW_LAYER_SIZE 0.05f
#define HAIR_SHADOW_FALL_OFF 5.0f
#define HAIR_SHADOW_BIAS 0.005f

layout(location = 0, index = 0) out vec4 color;

in vec3 hairDir;
in vec4 pos;

uniform mat4 objShadowMatrix;
uniform sampler2DShadow objShadow;

uniform mat4 hairShadowMatrix;
//uniform sampler2D hairShadow;
uniform sampler2D hairShadowDepth;

uniform vec3 lightDir;

void main()
{
	color.a = 1;
	vec3 dir = normalize(hairDir);

	//color = vec4(dir,1);
	//color.rgb = lightDir;
	//color.rgb = hairDir;
	//return;

	vec3 amb = vec3(0,0,0);
	vec3 aLight = AMBIENT_LIGHT;
	vec3 L = normalize(lightDir);
	float c = dot(dir, L);
	float s = sqrt( 1 - c*c );
	float dLight = LIGHT_INTENSITY * s;
	vec3 V = normalize(pos.xyz/pos.w);
	vec3 H = normalize(V + L);
	float sc = dot(dir, H);
	float ss = sqrt( 1 - sc*sc );
	float sLight = LIGHT_INTENSITY * s * pow( ss, SPECULAR_GLOSS );
	color.rgb = DIFFUSE_COLOR * dLight + sLight*SPECULAR_COLOR;
	amb = DIFFUSE_COLOR * aLight;

	vec4 objShadowPos = objShadowMatrix * pos;
	float visibility = texture(objShadow, vec3(objShadowPos.xy, objShadowPos.z/objShadowPos.w - SHADOW_BIAS));

	vec4 hairShadowPos = hairShadowMatrix * pos;
	float hairZ = texture(hairShadowDepth, hairShadowPos.xy).r;
	float z = hairShadowPos.z/hairShadowPos.w;
	float zDif = max( 0, ( z-HAIR_SHADOW_BIAS - hairZ ));

	visibility *= exp(-40.f*zDif);

	color.rgb *= visibility;

	color.rgb += amb;
	//color.rgb = vec3(0.68f,0.62f,0.18f);

	/*
	vec3 N = -normalize(hairDir);
    L = normalize(lightDir);
	vec3 ambient = vec3(0.2, 0.2, 0.2);
	vec3 diffuse = max(dot(N, L), 0.0) * vec3(0.68f,0.62f,0.18f) * 0.6;

	color = vec4(ambient + diffuse, 1.0);
	*/
}