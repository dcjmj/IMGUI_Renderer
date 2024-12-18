#version 330 core

layout(location = 0) in vec2 inVertex;

uniform mat4 cameraMatrix;
uniform sampler2D hairVerts;

void main()
{
	vec3 pos = texture2D(hairVerts,inVertex).xyz;
	//vec3 pos = texture2D(hairVerts,vec2(0,0)).xyz;
	gl_Position = cameraMatrix * vec4(pos,1);
	//gl_Position = vec4(inVertex,0,1);
}
