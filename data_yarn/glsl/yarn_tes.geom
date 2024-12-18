#version 430 core

layout (lines) in; 
 
//layout (line_strip, max_vertices = 256) out; 

#if 0

layout (points, max_vertices = 256) out;  
void main(void)         
{
	gl_Position = (gl_in[0].gl_Position);//+gl_in[1].gl_Position+gl_in[2].gl_Position)*0.33;
	gl_Position.xy = gl_Position.xy - vec2(0.5);
	EmitVertex();
	EndPrimitive();
}    

#else

//layout (points, max_vertices = 256) out; 
layout(triangle_strip, max_vertices=256) out;

#define SHADOW_BIAS 0.0001f
#define STEPS 18
#define M_PI 3.1415926535897932384626433832795

out float ran;
out vec3 yarnDir;
out vec4 pos;
out vec3 norm;
out vec4 shadow_coord;
out float distance2centralLine;
out vec2 texture_coord;
uniform mat4 shadow_matrix;
uniform mat4 cameraMatrix;
uniform float invCtrlPtsNum;
//uniform float segNum;
uniform sampler1D yarnVerts;
//uniform float rotations;

uniform sampler2D		color_tex;

float atan2(in float y, in float x)
{
    bool s = (abs(x) > abs(y));
    return mix(M_PI/2.0 - atan(x,y), atan(y,x), s);
}

/*
vec2 GetWDif( vec2 sPos1, vec2 sPos2 )
{
	vec2 sDir  = normalize(sPos2 - sPos1);
	vec2 sNorm = vec2( sDir.y, -sDir.x );
	return hair_width*sNorm;
}
*/

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

vec4 z =  vec4(0.0, 0.0, 1.0, 1.0);

void main()
{	
	// camera space
	vec4 cCenterPos1, cCenterPos2;
	vec4 width;
	vec3 dir;

	// model space
	vec4 mPos1, mPos2, mPos3, mPos4;
	vec4 mCenterPos1, mCenterPos2;

	vec2 cp1, cp2;
	vec3 d, aAxis, bAxis, offset, offset2;
	vec4 cMidPoint, cMidPoint2;
	float angle, radius, invSteps, rotations, s;

	float hair_width = 0.0025 * (12.0 / gl_in[0].gl_Position.w)* (12.0 / gl_in[0].gl_Position.w);

	vec4 samplePoint = gl_in[0].gl_Position;
	samplePoint.xy = samplePoint.xy - vec2(0.5);

	if (sqrt(samplePoint.x * samplePoint.x + samplePoint.y * samplePoint.y) < .5)
	{
		// get 4 ctrl pts
		float tc = (0.5 + int(gl_in[0].gl_Position.z)*3.0 )* invCtrlPtsNum;
		mPos1 = texture(yarnVerts, tc);
		tc += invCtrlPtsNum;
		mPos2 = texture(yarnVerts, tc);
		tc += invCtrlPtsNum;
		mPos3 = texture(yarnVerts, tc);
		tc += invCtrlPtsNum;
		mPos4 = texture(yarnVerts, tc);
	
		// configure angle and radius		
		radius = 0.08 * sqrt(samplePoint.x * samplePoint.x + samplePoint.y * samplePoint.y) * 2.0;
		invSteps = 1.0f/float(STEPS);
		rotations = texture(color_tex, vec2(samplePoint)).x * 6.0+10.0;	//float(STEPS);
		angle = atan2(samplePoint.x, samplePoint.y)+6.28/rotations*STEPS*int(gl_in[0].gl_Position.z);

		s = 0;

		mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
		cCenterPos2 = cameraMatrix * mCenterPos2;

		// BezierDeriv in camera space
		d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );

		aAxis = normalize( cross(d, z.xyz)  );
		bAxis = normalize( cross(aAxis, d) );

		offset2 = (aAxis * sin(angle) + bAxis * cos(angle)) * radius;

		cMidPoint2 = cCenterPos2 + vec4(offset2.xy, 0, 0);
		cMidPoint2.z += offset2.z * 0.01;

		// screen space
		cp2 = cCenterPos2.xy/cCenterPos2.z;

		// first n-1 points
		for (int step=0; step<STEPS; step++ ) {

			cMidPoint = cMidPoint2;
			cCenterPos1 = cCenterPos2;
			offset= offset2;

			s += invSteps;
			mCenterPos1 = mCenterPos2;
			mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
			cCenterPos2 = cameraMatrix * mCenterPos2;

			d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );

			aAxis = normalize( cross(d, z.xyz) );
			bAxis = normalize( cross(aAxis, d) );

			angle = angle+6.28/rotations;

			offset2 = normalize((aAxis * sin(angle)+ bAxis * cos(angle))) * radius;

			cMidPoint2 = cCenterPos2 + vec4( offset2.xy, 0, 0);
			cMidPoint2.z += offset2.z * 0.01;

			cp1 = cp2;
			cp2 = cCenterPos2.xy/cCenterPos2.z;

			dir = cCenterPos2.xyz + offset2 - cCenterPos1.xyz - offset;
			dir = normalize(dir);
			width = -vec4(normalize(cross(cMidPoint.xyz, dir).xy) * hair_width, 0.0, 0.0);
			gl_Position = cMidPoint - width;
			yarnDir = dir;
			pos = gl_Position;
			norm = offset;
			shadow_coord = shadow_matrix * mCenterPos1;
			shadow_coord.z -= SHADOW_BIAS;
			distance2centralLine = radius;
			texture_coord.x = s-invSteps+radius*10.0;
			texture_coord.y = 1.0;
			ran = radius*10.0;
			EmitVertex();
	
			gl_Position = cMidPoint + width;
			yarnDir = dir;
			pos = gl_Position;
			norm = offset;
			shadow_coord = shadow_matrix * mCenterPos1;
			shadow_coord.z -= SHADOW_BIAS;
			texture_coord.x = s-invSteps+radius*10.0;
			texture_coord.y = 0.0;
			distance2centralLine = radius;
			ran = radius*10.0;
			EmitVertex();
		}

		//------------- end points
		// get next 4 ctrl pts
	
		if (tc > 1.0 - 0.6*invCtrlPtsNum)	// prevent its coordinate out of range
		{
			mPos2 = mPos1; 
			mPos3 = mPos1; 
			mPos4 = mPos1; 
		}
		else
		{
			mPos1 = mPos4;
			tc += invCtrlPtsNum;
			mPos2 = texture(yarnVerts, tc);
			tc += invCtrlPtsNum;
			mPos3 = texture(yarnVerts, tc);
			tc += invCtrlPtsNum;
			mPos4 = texture(yarnVerts, tc);
		}

		s = 0;
		mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
		cCenterPos2 = cameraMatrix * mCenterPos2;

		d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );
		aAxis = normalize( cross(d, z.xyz)  );
		bAxis = normalize( cross(aAxis, d) );

		offset2 = (aAxis * sin(angle) + bAxis * cos(angle)) * radius;

		cMidPoint2 = cCenterPos2 + vec4(offset2.xy, 0, 0);
		cMidPoint2.z += offset2.z * 0.01;

		cMidPoint = cMidPoint2;
		cCenterPos1 = cCenterPos2;
		offset= offset2;

		mCenterPos1 = mCenterPos2;

		s = invSteps;
		mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
		cCenterPos2 = cameraMatrix * mCenterPos2;

		d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );
		aAxis = normalize( cross(d, z.xyz) );
		bAxis = normalize( cross(aAxis, d) );

		angle = angle+6.28/rotations;

		offset2 = (aAxis * sin(angle)+ bAxis * cos(angle)) * radius;

		cMidPoint2 = cCenterPos2 + vec4( offset2.xy, 0, 0);
		cMidPoint2.z += offset2.z * 0.01;

		cp1 = cMidPoint.xy/cMidPoint.z;
		cp2 = cMidPoint2.xy/cMidPoint2.z;

		dir = cCenterPos2.xyz + offset2 - cCenterPos1.xyz - offset;
		dir = normalize(dir);
		width = -vec4(normalize(cross(cMidPoint.xyz, dir).xy) * hair_width, 0.0, 0.0);

		gl_Position = cMidPoint - width;
		yarnDir = dir;
		pos = gl_Position;
		norm = offset;
		shadow_coord = shadow_matrix * mCenterPos1;
		shadow_coord.z -= SHADOW_BIAS;
		texture_coord.x = 1+radius*10.0;
		texture_coord.y = 1.0;
		distance2centralLine = radius;
		ran = radius*10.0;
		EmitVertex();
	
		gl_Position = cMidPoint + width;
		yarnDir = dir;
		pos = gl_Position;
		norm = offset;
		shadow_coord = shadow_matrix * mCenterPos1;
		shadow_coord.z -= SHADOW_BIAS;
		texture_coord.x = 1+radius*10.0;
		texture_coord.y = 0.0;
		distance2centralLine = radius;
		ran = radius*10.0;
		EmitVertex();
	
		EndPrimitive();
	}
}
#endif