#version 420 core

layout(location = 0) in vec3 inVertex;

uniform mat4 cameraMatrix;
uniform sampler1D ctlPts;
uniform float ctlPtsNum;

//out vec3 vPosition;

void main()
{
	float tc = (.5+inVertex.x)/ctlPtsNum;
	vec3 pos = texture(ctlPts, tc).xyz;
	gl_Position = cameraMatrix * vec4(pos,1);
	//gl_Position = vec4(inVertex,1);

	//vPosition = inVertex;
	
	//const vec4 vertices[] = vec4[](vec4( 0.4, -0.4, 0.5, 1.0), vec4(-0.4, -0.4, 0.5, 1.0), vec4( 0.4,  0.4, 0.5, 1.0), vec4(-0.4,  0.4, 0.5, 1.0));                                                              
    //gl_Position = vertices[gl_VertexID];
}