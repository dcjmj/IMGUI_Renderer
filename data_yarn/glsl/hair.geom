#version 330 core

#define HAIR_WIDTH 0.01f
#define STEPS 6

layout(points) in;
layout(triangle_strip, max_vertices=93) out;

out vec3 hairDir;
out vec4 pos;

uniform mat4 cameraMatrix;
uniform mat3 cameraNormalMatrix;
uniform float xskip;
uniform sampler2D hairVerts;

vec2 GetWDif( vec2 sPos1, vec2 sPos2 )
{
	vec2 sDir  = normalize(sPos2 - sPos1);
	vec2 sNorm = vec2( sDir.y, -sDir.x );
	return HAIR_WIDTH*sNorm;
}

vec4 CatmullRomTangentDir( vec4 cPos0, float s0, vec4 cPos1, float s1, vec4 cPos2, float s2 )
{
	float d01 = s1 - s0;
	float d12 = s2 - s1;
	return ( (d12*d12)*(cPos1-cPos0) + (d01*d01)*(cPos2-cPos1) ) / (d01+d12) * (1.0f/3.0f);
}

vec3 CatmullRomTangentDir( vec3 pos0, float s0, vec3 pos1, float s1, vec3 pos2, float s2 )
{
	float d01 = s1 - s0;
	float d12 = s2 - s1;
	return ( (d12*d12)*(pos1-pos0) + (d01*d01)*(pos2-pos1) ) / (d01+d12) * (1.0f/3.0f);
}

vec4 Bezier( float s, vec4 b0, vec4 b1, vec4 b2, vec4 b3 )
{
	float t = 1-s;
	float tt = t*t;
	float ts = t*s;
	float ss = s*s;
	return tt*t*b0 + 3*tt*s*b1 + 3*t*ss*b2 + s*ss*b3;
}

vec3 BezierDeriv( float s, vec3 b0, vec3 b1, vec3 b2, vec3 b3 )
{
	float t = 1-s;
	return 3*t*t*(b1-b0) + 6*t*s*(b2-b1) + 3*s*s*(b3-b2);
}

void main()
{
	vec2 tc = gl_in[0].gl_Position.xy;

	vec3 pos1 = texture2D(hairVerts,tc).xyz;
	vec4 cPos1 = cameraMatrix * vec4(pos1,1);
	vec2 sPos1 = cPos1.xy / cPos1.z;
	vec4 bezier1 = cPos1;
	vec3 bez1 = pos1;

	vec4 cPos0 = cPos1;
	vec4 bezier0 = cPos0;
	vec3 bez0 = bez1;

	tc.x += xskip;
	vec3 pos2 = texture2D(hairVerts,tc).xyz;
	vec4 cPos2 = cameraMatrix * vec4(pos2,1);
	vec2 sPos2 = cPos2.xy / cPos2.z;

	hairDir = vec3(0,0,0);

	const float s0=0, s1=1, s2=2, s3=3;	// uniform parameterization

	tc.x += xskip;
	while ( tc.x < 1 ) {
		vec3 pos3 = texture2D(hairVerts,tc).xyz;
		vec4 cPos3 = cameraMatrix * vec4(pos3,1);	// camera space
		vec2 sPos3 = cPos3.xy / cPos3.z;

		vec4 tangentDir = CatmullRomTangentDir( cPos1, s1, cPos2, s2, cPos3, s3 );
		vec4 tangentIn  = tangentDir / (s2-s3);
		vec4 tangentOut = tangentDir / (s2-s1);

		vec4 bezier2 = cPos2 + tangentIn;
		vec4 bezier3 = cPos2;

		vec3 ptangentDir = CatmullRomTangentDir( pos1, s1, pos2, s2, pos3, s3 );
		vec3 ptangentIn  = ptangentDir / (s2-s3);
		vec3 ptangentOut = ptangentDir / (s2-s1);

		vec3 bez2 = pos2 + ptangentIn;
		vec3 bez3 = pos2;

		vec2 wdif = GetWDif(sPos1, sPos2);	// screen space

		gl_Position = cPos1;
		pos = gl_Position;
		EmitVertex();
		gl_Position.xy += wdif;
		pos = gl_Position;
		EmitVertex();

		for ( int step=1; step<STEPS; step++ ) {
			float s = float(step)/float(STEPS);
			gl_Position = Bezier( s, bezier0, bezier1, bezier2, bezier3 );
			hairDir = cameraNormalMatrix * normalize(BezierDeriv( s, bez0, bez1, bez2, bez3 ));
			pos = gl_Position;
			EmitVertex();
			gl_Position.xy += wdif;
			pos = gl_Position;
			EmitVertex();
		}

		tc.x += xskip;
		cPos0 = cPos1;
		cPos1 = cPos2;
		cPos2 = cPos3;
		sPos1 = sPos2;
		sPos2 = sPos3;

		bezier0 = bezier3;
		bezier1 = bezier3 + tangentOut;
		bez0 = bez3;
		bez1 = bez3 + ptangentOut;
		hairDir = cameraNormalMatrix * normalize(ptangentOut);
		pos1 = pos2;
		pos2 = pos3;
	}

	vec2 wdif = GetWDif(sPos1, sPos2);

	gl_Position = cPos1;
	pos = gl_Position;
	EmitVertex();
	gl_Position.xy += wdif;
	pos = gl_Position;
	EmitVertex();

	vec4 bezier2 = cPos2;
	vec4 bezier3 = cPos2;
	vec3 bez2 = pos2;
	vec3 bez3 = pos2;
	for ( int step=1; step<=STEPS; step++ ) {
		float s = float(step)/float(STEPS);
		gl_Position = Bezier( s, bezier0, bezier1, bezier2, bezier3 );
		hairDir = cameraNormalMatrix * BezierDeriv( s, bez0, bez1, bez2, bez3 );
		pos = gl_Position;
		EmitVertex();
		gl_Position.xy += wdif;
		pos = gl_Position;
		EmitVertex();
	}

	EndPrimitive();

	//gl_InvocationID
}