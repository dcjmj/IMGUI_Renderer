#version 330 core

//layout(location = 0, index = 0) out vec4 color;
layout(location = 0, index = 0) out vec4 depth;

in vec4 pos;

void main()
{
	//color = vec4(gl_FragCoord.z);
	//color = vec4(0, 0, 1, 1);

	depth = vec4(gl_FragCoord.z, gl_FragCoord.z, 1, 1);
}