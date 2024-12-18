#version 330 core

layout(location = 0, index = 0) out vec4 color;

vec3 defaultColor = vec3(0.2f, 0.8f,0.1f);

in vec3 norm;
in vec3 pos;
uniform vec3 light_pos_world;
uniform mat4 camera_matrix;

void main()
{
	color = vec4(1, 0, 0, 1);

	//vec3 N = vec3(norm.x, norm.y, norm.z); 
	//vec3 L = -normalize(light_pos_world + pos);

	//color.xyz = max(dot(L, normalize(N) ), 0) * defaultColor;

	//color.xyz =  normalize(pos);
}