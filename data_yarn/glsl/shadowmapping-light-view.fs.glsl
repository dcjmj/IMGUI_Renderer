#version 420

uniform sampler2D tex_depth;

layout (location = 0) out vec4 color;

void main(void)
{
    float d = texelFetch(tex_depth, ivec2(gl_FragCoord.xy * 3.0) + ivec2(850, 1050), 0).r;
    d = (d - 0.95) * 5.0;
    color = vec4(d, d, d, d);
	//color = vec4(1.0,0.0,0.0,1.0);
}
