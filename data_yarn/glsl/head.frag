#version 460 core

layout(location = 0) out vec4 color;

in vec3 vertex_normal;
in vec3 posWorld;
uniform sampler2DShadow shadow_tex;
uniform mat4 shadow_matrix;
uniform vec3 cameraPos;
uniform vec3 light_dir;
#define SHADOW_BIAS 0.002f      

vec3 gammaCorrection (vec3 colour, float gamma) {
  return pow(colour, vec3(1. / gamma));
}

float PCF(vec4 shadow_tex_coords, sampler2DShadow shadow_tex, int kernel_size, float texel_size) {
    float shadow = 0.0;
    int half_kernel = kernel_size / 2;
    for (int x = -half_kernel; x <= half_kernel; ++x) {
        for (int y = -half_kernel; y <= half_kernel; ++y) {
            vec4 offset_coords = shadow_tex_coords;
            offset_coords.xy += vec2(x, y) * shadow_tex_coords.w * texel_size;
            shadow += textureProj(shadow_tex, offset_coords).r;
        }
    }
    return shadow / (kernel_size * kernel_size);
}

void main()
{
	bool blinn = false;
	vec3 color_based = vec3(0.25,0.24,0.22);
    // ambient
    vec3 ambient = 0.05 * color_based;
    // diffuse
	vec3 color_ = vec3(0,0,0);
 	{
		vec3 FragColor =vec3(0,0,0);

        vec3 L = vec3(4, 4, 4);
        vec3 lightDir = normalize(light_dir);
        bool direct = true;
        vec3 normal = normalize(vertex_normal);
        //diffuse
        float diff = max(dot(lightDir, normal), 0.0);
        vec3 diffuse = diff * color_based;
        // specular
        vec3 viewDir = normalize(cameraPos - posWorld);
        float spec = 0.0;
        if(blinn)
        {
            vec3 halfwayDir = normalize(lightDir + viewDir);  
            spec = pow(max(dot(normal, halfwayDir), 0.0), 32.0);
        }
        else
        {
            vec3 reflectDir = reflect(-lightDir, normal);
            spec = pow(max(dot(viewDir, reflectDir), 0.0), 8.0);
        }
        vec3 specular = vec3(0.07) * spec; // assuming bright white light color
        specular = vec3(0);
        //vec3 specular = vec3(0.f,0.f,0.f);
        vec3 contrib = (ambient + diffuse + specular) * L;
        FragColor = FragColor + contrib;

		float shadow = 0.0;
		
        float texel_size = 1.0 / textureSize(shadow_tex, 0).x;  // 假设方形纹理
        int kernel_size = 3;  // 3x3 的核
        vec4 real_world_pos = vec4(posWorld,1);
        vec4 world_pos_here;
        world_pos_here	= shadow_matrix * real_world_pos;

        world_pos_here.z -= SHADOW_BIAS;

        float global_shadow_scale = PCF(world_pos_here, shadow_tex, kernel_size, texel_size);
        
		FragColor = FragColor * global_shadow_scale;
        
		color_+=FragColor;
		//shadow = 0.0;
	}
	//color_ = gammaCorrection(color_, 2.2);
	color = vec4(color_, 1.0);
};