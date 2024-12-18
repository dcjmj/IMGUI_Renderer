#version 420

uniform sampler2D depthbuffer_tex;

in vec2 UV;

layout (location = 0) out vec4 color;

void main(void)
{
	color = vec4( texture( depthbuffer_tex, UV).w );
	color = vec4( texture( depthbuffer_tex, UV).r );

	color = texture( depthbuffer_tex, UV);
	//color.xyz = vec3(UV, 0.0f);

	//float d = texelFetch(depthbuffer_tex, ivec2(gl_FragCoord.xy * 1.0) + ivec2(1024, 1024), 0).r;
    //d = (d - 0.99)*30.0;
	color = vec4( ( texture( depthbuffer_tex, UV).r- 0.995)*200.0 );
    //color = vec4(d);
}
