#version 330 core

layout(location = 0, index = 0) out vec4 color;

#define HAIR_WIDTH 0.0025f

in float ran;
in vec3 yarnDir;
in vec4 pos;
in vec3 norm;
in vec4 shadow_coord;
in float distance2centralLine;
in vec2 texture_coord;

uniform bool enable_shadow;

uniform sampler2DShadow shadow_tex;
//uniform sampler2D		fiber_norm_tex;
//uniform sampler2D		color_tex;

uniform vec3 lightPos;
uniform vec3 specular_albedo = vec3(0.26);
uniform float specular_power = 128.0;

uniform mat4 cameraMatrix;

float specular, diffuse, ambient;
//vec3 ambient = vec3(0.02, 0.08, 0.01);
vec3 V, N, L;

vec4 defaultColor = vec4(0.2f,0.8f,0.1f, 1.0f);
//uniform bool colorful;
void main()
{
	// texture
	vec2 tex_c = texture_coord;
	tex_c.x *= 4.0;

	vec3 local_yarn_dir = yarnDir;

	ambient = 0.1;
	color.a = 1.0;
	//vec3 dir = normalize(local_yarn_dir);
	//color = vec4(dir,1);
	//color.rgb = lightPos;
	//color.rgb = normalize(local_yarn_dir);
	//color.rgb = vec3(0.68f,0.62f,0.18f);

	// diffuse
	N = normalize(norm);
	
    L = normalize(lightPos - pos.xyz);

	float c = dot(local_yarn_dir, L);
	float s = sqrt( 1 - c*c );

	diffuse = (s * 0.75 + 0.25);

	// specular
	V = -normalize(pos.xyz/pos.w);
	vec3 H = normalize(V + L);
	float sc = dot(local_yarn_dir, H);
	float ss = sqrt( 1 - sc*sc );

	specular = s * pow(ss, specular_power )  ;

	// height mapping
	
	//float Hgiven = texture(fiber_norm_tex, tex_c).x;
	//tex_c.x += 0.02;
	//float Hright = texture(fiber_norm_tex, tex_c).x; 
	//tex_c.x -= 0.02;
	//tex_c.y += 0.02;
	//float Habove = texture(fiber_norm_tex, tex_c).x;
	
	//vec3 tmpN = -normalize(vec3(Hgiven-Habove, Hgiven-Hright, 1.0));

	// visibility
	// max(dot(-texture(fiber_tex, tex_c).xyz, L), .0) *
	//float visibility = max(dot(N-0.5*HAIR_WIDTH*texture(fiber_norm_tex, tex_c).xyz, L), .0)  ;//* max(dot(N, L), .0) * 0.6;
	//float visibility = max(dot(normalize( mat3(cameraMatrix)*texture(fiber_norm_tex, tex_c).xyz ), L), .0) ;//* max(dot(N, L), 0.0)  ;
	//float visibility = max(dot(tmpN, L), 0.0) ;//+ max(dot(N, L), 0.0)*0.2;

	float visibility = max(dot(N, L), 0.0);
	vec3 z = vec3(0,0,1);
	vec3 zxd = normalize(cross(z, yarnDir));
	N = -normalize(sqrt(1-(0.5-tex_c.y)*(0.5-tex_c.y)) * z + (0.5-tex_c.y) * zxd);
	diffuse *= max(dot(N, L), 0.0);

	// shadow
	float inShadow = enable_shadow ? textureProj(shadow_tex, shadow_coord) : 1.0;
	float ambientOcclusion = max(distance2centralLine * 10.0f - 0.15f, 0.0f) * 2.5f;
	//float tmp = ambient + inShadow * (visibility) * (diffuse + specular)* ambientOcclusion ;//* (diffuse + specular);
	float tmp =   inShadow * 0.8 * ambientOcclusion * (diffuse ) ;// visibility *;

	// texture
	
	//if (colorful)
	//	defaultColor.xyz = texture(color_tex, vec2(ran, ran)).xyz * texture(fiber_norm_tex, tex_c).xyz;
	//else
	//	defaultColor.xyz *= texture(fiber_norm_tex, tex_c).xyz;
	//if (defaultColor.x < 0.1) discard;
	
	color = tmp * defaultColor;
	color.w = 1.0;

	//color.xyz = mat3(cameraMatrix)*texture(fiber_norm_tex, tex_c).xyz ;
	
	//color.xyz = tmpN;

	//color.xyz = texture(fiber_norm_tex,tex_c).xyz;

	//defaultColor = texture(fiber_tex,tex_c);
	
	/*
	if (defaultColor.x > 0.2) {
		//defaultColor.a = 0.0;
		
		color = tmp * defaultColor;
		color.w = 1.0;
	} else {
		color.rgb = vec3(1.0, 0.0, 0.0);
		color.w = 0.0;
	}
	*/

	//color = tmp * defaultColor;
	//color.a = 1.0;

	//color.rgb = vec3(inShadow);

	//vec3 dir = normalize(local_yarn_dir);
	//color = vec4(dir,1);
	//color = vec4(normalize(norm),1);
	//color = vec4(lightPos,1);
	//color.rgb = max(dot(N, L), 0.0) * vec3(0.2f,0.8f,0.1f) * 0.6;
	//color = vec4(gl_FragCoord.z);

	//color = vec4(texture_coord, 0.0, 1.0);
	//color = texture(fiber_tex,texture_coord);

	//color.xyz = vec3(sin(tex_c.y),cos(tex_c.y),1);
}