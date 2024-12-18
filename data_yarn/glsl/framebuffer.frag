#version 420

uniform sampler2D framebuffer_tex;

in vec2 UV;

layout (location = 0) out vec4 color;

void main(void)
{
    //float d = texelFetch(framebuffer_tex, ivec2(gl_FragCoord.xy * 3.0) + ivec2(850, 1050), 0).r;
    //d = (d - 0.99)*30.0;
    //color = vec4(d);

	//color = vec4(texture( framebuffer_tex, (vec2(1.0) + gl_FragCoord.xy) * 0.5f ).xyz, 1.0);
	
	//color = vec4(texelFetch(framebuffer_tex, ivec2(gl_FragCoord.xy * 1.0) + ivec2(1280, 960), 0).xyz, 1.0);
	
	//color = vec4(UV, 0.0, 1.0);
	color = vec4(texture( framebuffer_tex, UV).xyz * 2.0f - vec3(1.0, 1.0, 1.0), 1.0);
	//color = vec4(1.0, 0.0, 0.0, 1.0);
	//color = vec4(texture( framebuffer_tex, UV).w) * 2.0 - 1.0;
}
