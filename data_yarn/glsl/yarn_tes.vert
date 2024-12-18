#version 430 core

uniform float segmentNum;
out float vSegmentNum;

void main()
{                                                     
    const vec4 vertices[] = vec4[](vec4( 0.1, -0.1, 0.5, 1.0),    
                                   vec4(-0.1, -0.1, 0.5, 1.0),    
                                   vec4(-0.1,  0.1, 0.5, 1.0),    
                                   vec4( 0.1,  0.1, 0.5, 1.0)); 
	                                                      
    gl_Position = vertices[gl_VertexID];   
	vSegmentNum = segmentNum;
}
