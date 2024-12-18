#version 330 core

layout (location = 0) in vec3 position;   // 输入的顶点坐标

out defaultBlock
{
	vec4 position;
	vec2 uv;
} outBlock;

void main()
{
	outBlock.position = vec4(position, 1.0);
	outBlock.uv = outBlock.position.xy * 0.5f + 0.5f;
	gl_Position = outBlock.position;
}