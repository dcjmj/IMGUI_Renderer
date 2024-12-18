#version 330 core

#define HAIR_WIDTH 0.0025f
#define SHADOW_BIAS 0.0001f
#define STEPS 17

layout(points) in;

layout(triangle_strip, max_vertices=68) out;
//layout(points, max_vertices=68) out;
//layout(line_strip, max_vertices=68) out;

int segmentCount = 8;
float expasion = 0.05;

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
uniform float segNum;
uniform sampler1D yarnVerts;
//uniform float rotations;

vec2 GetWDif( vec2 sPos1, vec2 sPos2 )
{
	vec2 sDir  = normalize(sPos2 - sPos1);
	vec2 sNorm = vec2( sDir.y, -sDir.x );
	return HAIR_WIDTH*sNorm;
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

vec4 z =  vec4(0.0, 0.0, 1.0, 0.0);

void main()
{
	// camera space
	vec4 cCenterPos1, cCenterPos2;
	vec4 width;
	vec3 dir;

	// model space
	vec4 mPos1, mPos2, mPos3, mPos4;
	vec4 mCenterPos1, mCenterPos2;

	// get 4 ctrl pts
	float tc = (0.5 + int(gl_in[0].gl_Position.w*segNum)*3.0 )* invCtrlPtsNum;
	mPos1 = texture(yarnVerts, tc);

	tc += invCtrlPtsNum;
	mPos2 = texture(yarnVerts, tc);

	tc += invCtrlPtsNum;
	mPos3 = texture(yarnVerts, tc);

	tc += invCtrlPtsNum;
	mPos4 = texture(yarnVerts, tc);

	float invSteps = 1.0f/float(STEPS);
	float s = (gl_in[0].gl_Position.w*segNum - int(gl_in[0].gl_Position.w*segNum))*4.0;
	float angle = gl_in[0].gl_Position.z;

	/*
	vec4 tmp1 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
	vec4 tmp2 = Bezier(s+invSteps, mPos1, mPos2, mPos3, mPos4);
	vec4 direction = normalize(tmp2 - tmp1);
	
	mPos1 = tmp1;
	mPos2 = mPos1 - z * 0.01 * log(abs(cos(angle)));
	mPos3 = mPos2 - z * 0.01 * log(abs(cos(angle)));
	mPos4 = mPos3 - z * 0.01 * log(abs(cos(angle)));
	*/
	//z =  vec4(1.0, 0.0, 0.0, 1.0);

	// configure angle and radius
	//float rotations = 18;
	
	float radius = gl_in[0].gl_Position.x;
	
	float rotations =  gl_in[0].gl_Position.y;
	//float s = 0;
	vec3 d, aAxis, bAxis, offset, offset2;
	vec4 cMidPoint, cMidPoint2;

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
	vec2 cp1;
	vec2 cp2 = cCenterPos2.xy/cCenterPos2.z;
	//s = 0;
	// first n-1 points
	float r;
	int count = 0;
	for (int step=0; step<segmentCount; step++ ) {
		
		if (s > 1.0-invSteps)
			break;
		count++;
		r = radius +  (radius)*expasion*count;

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

		offset2 = normalize((aAxis * sin(angle)+ bAxis * cos(angle))) * r;

		cMidPoint2 = cCenterPos2 + vec4( offset2.xy, 0, 0);
		cMidPoint2.z += offset2.z * 0.01;

		cp1 = cp2;
		cp2 = cCenterPos2.xy/cCenterPos2.z;

		dir = cCenterPos2.xyz + offset2 - cCenterPos1.xyz - offset;
		dir = normalize(dir);
		width = -vec4(normalize(cross(cMidPoint.xyz, dir).xy) * HAIR_WIDTH, 0.0, 0.0);
		gl_Position = cMidPoint + width;
		yarnDir = dir;
		pos = gl_Position;
		norm = offset;
		shadow_coord = shadow_matrix * mCenterPos1;
		shadow_coord.z -= SHADOW_BIAS;
		distance2centralLine = r;
		texture_coord.x = s-invSteps;
		texture_coord.y = 0;
		ran = gl_in[0].gl_Position.x*10.0;
		EmitVertex();
	
		gl_Position = cMidPoint - width;
		yarnDir = dir;
		pos = gl_Position;
		norm = offset;
		shadow_coord = shadow_matrix * mCenterPos1;
		shadow_coord.z -= SHADOW_BIAS;
		texture_coord.x = s-invSteps;
		texture_coord.y = 1.0;
		distance2centralLine = r;
		ran = gl_in[0].gl_Position.x*10.0;
		EmitVertex();
	}

	//------------- end points
	// get next 4 ctrl pts
	
	if (count < segmentCount) {
		s = s - 1.0 + invSteps;

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
		mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
		cCenterPos2 = cameraMatrix * mCenterPos2;

		// BezierDeriv in camera space
		d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );

		aAxis = normalize( cross(d, z.xyz)  );
		bAxis = normalize( cross(aAxis, d) );

		offset2 = (aAxis * sin(angle) + bAxis * cos(angle)) * (radius + (radius)*expasion*count);

		cMidPoint2 = cCenterPos2 + vec4(offset2.xy, 0, 0);
		cMidPoint2.z += offset2.z * 0.01;

		// screen space
		vec2 cp1;
		vec2 cp2 = cCenterPos2.xy/cCenterPos2.z;

		// first n-1 points
		float r;
		for (int step=0; step<segmentCount; step++ ) {
			if (count > segmentCount)
				break;
			count++;

			r = radius + (radius)*expasion*count;

			cMidPoint = cMidPoint2;
			cCenterPos1 = cCenterPos2;
			offset= offset2;

			s += invSteps;

			if (s > 1.0)
				break;

			mCenterPos1 = mCenterPos2;
			mCenterPos2 = Bezier(s, mPos1, mPos2, mPos3, mPos4);
			cCenterPos2 = cameraMatrix * mCenterPos2;

			d = normalize( mat3( cameraMatrix) * BezierDeriv(s, mPos1.xyz, mPos2.xyz, mPos3.xyz, mPos4.xyz) );

			aAxis = normalize( cross(d, z.xyz) );
			bAxis = normalize( cross(aAxis, d) );

			angle = angle+6.28/rotations;

			offset2 = normalize((aAxis * sin(angle)+ bAxis * cos(angle))) * r;

			cMidPoint2 = cCenterPos2 + vec4( offset2.xy, 0, 0);
			cMidPoint2.z += offset2.z * 0.01;

			cp1 = cp2;
			cp2 = cCenterPos2.xy/cCenterPos2.z;

			dir = cCenterPos2.xyz + offset2 - cCenterPos1.xyz - offset;
			dir = normalize(dir);
			width = -vec4(normalize(cross(cMidPoint.xyz, dir).xy) * HAIR_WIDTH, 0.0, 0.0);
			gl_Position = cMidPoint + width;
			yarnDir = dir;
			pos = gl_Position;
			norm = offset;
			shadow_coord = shadow_matrix * mCenterPos1;
			shadow_coord.z -= SHADOW_BIAS;
			distance2centralLine = r;
			texture_coord.x = s+1.0-invSteps;
			texture_coord.y = 0;
			ran = gl_in[0].gl_Position.x*10.0;
			EmitVertex();
	
			gl_Position = cMidPoint - width;
			yarnDir = dir;
			pos = gl_Position;
			norm = offset;
			shadow_coord = shadow_matrix * mCenterPos1;
			shadow_coord.z -= SHADOW_BIAS;
			texture_coord.x = s+1.0-invSteps;
			texture_coord.y = 1.0;
			distance2centralLine = r;
			ran = gl_in[0].gl_Position.x*10.0;
			EmitVertex();
		}
	}
	
	EndPrimitive();
}