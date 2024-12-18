#version 420 core                                              
        
out vec2 UV;
		                                                       
void main(void)                                                
{                                                              
    const vec4 vertices[] = vec4[](vec4(-1.0, -1.0, 0.5, 1.0), 
                                   vec4( 1.0, -1.0, 0.5, 1.0), 
                                   vec4(-1.0,  1.0, 0.5, 1.0), 
                                   vec4( 1.0,  1.0, 0.5, 1.0));
                                                               
    gl_Position = vertices[gl_VertexID];
	UV = (vec2(1.0) + vertices[gl_VertexID].xy) * 0.5f;                  
}       
