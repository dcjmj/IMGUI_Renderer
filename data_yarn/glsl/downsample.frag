#version 330 core

in defaultBlock
{
	vec4 position;
	vec2 uv;
} inBlock;

out vec4 FragColor;  // 输出颜色

uniform sampler2D uTexture;  // 输入的 2Ax2B 的纹理
uniform sampler2D uTexture2;  // 输入的 2Ax2B 的纹理

void main() {
    vec2 texSize = textureSize(uTexture, 0);  // 返回纹理的宽高

    vec2 texelSize = vec2(1.0) / texSize;     // 每个纹理像素的大小

    // find the topleft pixel
    vec2 topLeft = inBlock.uv;

    // get the color of each
    vec4 color1 = texture(uTexture, topLeft);
    vec4 color2 = texture(uTexture, topLeft + vec2(texelSize.x, 0.0));
    vec4 color3 = texture(uTexture, topLeft + vec2(0.0, texelSize.y));
    vec4 color4 = texture(uTexture, topLeft + texelSize);
   // vec4 color5 = texture(uTexture, topLeft + texelSize * 2);
    //vec4 color6 = texture(uTexture, topLeft + vec2(texelSize.x * 2, 0.0));
   // vec4 color7 = texture(uTexture, topLeft + vec2(0.0, texelSize.y * 2));
    //vec4 color8 = texture(uTexture, topLeft + vec2(texelSize.x, texelSize.y * 2));
   // vec4 color9 = texture(uTexture, topLeft + vec2(texelSize.x * 2, texelSize.y));

    //FragColor = (color1 + color2 + color3 + color4 + color5 +color6 +color7 +color8+color9) / 9.0;
    FragColor = (color1 + color2 + color3 + color4) / 4.0;
    FragColor.w = 1.0;
}
