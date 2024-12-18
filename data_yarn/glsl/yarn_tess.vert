#version 430 core

layout(location = 0) in vec4 inVertex;

out float arclength;
out int visible;

uniform mat4 view_matrix;

bool isVisible(vec3 vert)
{
    vec4 p = view_matrix * vec4(vert,1);
    return !(( p.x < -p.w)||
			 ( p.x >  p.w)||
             ( p.y < -p.w)||
             ( p.y >  p.w)||
             ( p.z < 0));
}

void main()
{    
	gl_Position = vec4(inVertex.xyz, 1.0f);
	arclength = inVertex.w;

	visible =  isVisible(inVertex.xyz) ? 1 : 0;
}
