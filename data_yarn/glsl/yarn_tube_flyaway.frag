#version 430 core

#define M_PI		3.141592653589793238462643383279f
#define M_PI_2		6.283185307179586476925286766559f
#define M_PI_2_INV	0.159154943091895335768883763372f
#define SHADOW_BIAS 0.0001f

uniform vec3 AR = vec3(0.15, 0.15, 0.15);
uniform vec3 ATT = vec3(0.15, 0.3, 0.21);
uniform vec3 AD = vec3(0.09, 0.35, 0.135);
uniform float betaR = 0.6;
uniform float betaTT = 0.2;
uniform float betaN = 0.6;
uniform float d_b = 0.8;
uniform float d_f = 0.8;
uniform float beta_diffuse = 0.902;
uniform float d_f_inner = 0.8;
uniform float d_cross_inter = 0.2;

float atan2(in float y, in float x)
{
    bool s = (abs(x) > abs(y));
    return mix(M_PI/2.0 - atan(x,y), atan(y,x), s);
}

float erf(float x)
{
    // constants
    float a1 =  0.254829592;
    float a2 = -0.284496736;
    float a3 =  1.421413741;
    float a4 = -1.453152027;
    float a5 =  1.061405429;
    float p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = abs(x);

    // A&S formula 7.1.26
    float t = 1.0/(1.0 + p*x);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

struct Frame{
	vec3 t, b, n;
};

vec3 to_local(Frame frame, vec3 w) { 
	return vec3(dot(frame.t, w), dot(frame.b, w), dot(frame.n, w)); 
}

vec3 to_world(Frame frame, vec3 w) 
{
	return frame.t * w.x + frame.b * w.y + frame.n * w.z; 
}
float trigInverse(float x) { return min(sqrt(max(1.0f - x * x, 0.0f)), 1.0f); }

float integrateG(float x1, float x2, float sigma)
{
	return 0.5f * (erf(x2 / (sqrt(2.0f) * sigma)) - erf(x1 / (sqrt(2.0f) * sigma)));
}

float integrateG_phi(float sigma)
{
	float result = 0.f;
	float delta;
	float shift = 0.0f;
	do {
		delta = integrateG(-M_PI / 2.0f + shift, M_PI / 2.0f + shift, sigma) +
			integrateG(-M_PI / 2.0f - shift - 2 * M_PI, M_PI / 2.0f - shift - 2 * M_PI, sigma);
		result += delta;
		shift += 2.0 * M_PI;
	} while (delta > 1e-4f);
	return result;
}
vec3 fR(float theta_i) 
{
	return AR + (vec3(1.0f, 1.0f, 1.0f) - AR) * pow(1.0f - cos(theta_i), 5.0f);
}

vec3 getAlphaForward_Single(float theta_i, float phi_i, int lobe) 
{
	float M, N;
	float cosThetaI = clamp(cos(theta_i), 0.f, 1.f);
	vec3 a;
	switch (lobe) {
	case 0:
		// M = integrateG(-M_PI / 2.0f + theta_i, M_PI / 2.0f + theta_i, _betaR);
		M = 1.0f;
		N = 0.5f;
		// a = M * N * fR(theta_i) * cos(theta_i);
		a = M * N * fR(theta_i);
		break;
	case 1:
		M = 1.0f;
		N = integrateG_phi(betaN);
		a = M * N * ATT * (vec3(1, 1, 1) - fR(theta_i));
		break;
	case 2:
		M = 1.0f;
		N = 0.5f;
		// a = M * N * AD * cos(theta_i);
		a = M * N * AD;
		break;
	}
	return a * cosThetaI;
}

vec3 getAlphaBackward_Single(float theta_i, float phi_i, int lobe)
{
	float M, N;
	vec3 a;
	float cosThetaI = clamp(cos(theta_i), 0.f, 1.f);
	switch (lobe) {
	case 0:
		// M = integrateG(-M_PI / 2.0f + theta_i, M_PI / 2.0f + theta_i, _betaR);
		M = 1.0f;
		N = 0.5f;
		// a = M * N * fR(theta_i) * cos(theta_i);
		a = M * N * fR(theta_i);
		break;
	case 1:
		// M = integrateG(-M_PI / 2.0f + theta_i, M_PI / 2.0f + theta_i, _betaTT);
		M = 1.0f;
		N = integrateG_phi(betaN);
		N = 1 - N;
		// a = M * N * ATT * (vec3(1, 1, 1) - fR(theta_i)) * cos(theta_i);
		a = M * N * ATT * (vec3(1, 1, 1) - fR(theta_i));
		break;
	case 2:
		M = 1.0f;
		N = 0.5f;
		// a = M * N * AD * cos(theta_i);
		a = M * N * AD;
		break;
	}
	return a * cosThetaI;
}

float g_detector(float beta, float theta)
{
	return exp(-theta * theta / (2.0f * beta * beta)) / (sqrt(2.0f * M_PI) * beta);
}

float g_theta(float beta, float phi)
{
	float result = 0.0f;
	float delta;
	float shift = 0.0f;
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 1 * M_PI);
		result += delta;
		shift += 1 * M_PI;
	}
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 1 * M_PI);
		result += delta;
		shift += 1 * M_PI;
	}
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 1 * M_PI);
		result += delta;
		shift += 1 * M_PI;
	}
	return result;
}

float g_phi(float beta, float phi)
{
	float result = 0.0f;
	float delta;
	float shift = 0.0f;
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 2 * M_PI);
		result += delta;
		shift += 2 * M_PI;
	}
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 2 * M_PI);
		result += delta;
		shift += 2 * M_PI;
	}
	{
		delta = g_detector(beta, phi + shift) + g_detector(beta, phi - shift - 2 * M_PI);
		result += delta;
		shift += 2 * M_PI;
	}
	return result;
}

vec3 getAlphaForward(float theta_i, float phi_i, int lobe, vec3 v_F_M_ply, vec3 v_F_N_ply,
vec3 v_B_M_ply, vec3 a_F_ply, vec3 a_B_ply)
{
	float M;
	vec3 N;
	float cosThetaI = clamp(cos(theta_i), 0.f, 1.f);
	vec3 a;
	float theta_I_fiber = asin(sin(theta_i) * 1.f);
	switch (lobe) {
	case 0:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, fiber_bsdf->_betaR);
		M = 1.f;
		N = vec3(0.f);
		a = M * N * fR(theta_I_fiber);
		// a = a * cos(theta_I_fiber);
		break;
	case 1:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, _betaTT);
		M = 1.f;
		N = vec3(integrateG_phi(v_F_N_ply[0]), integrateG_phi(v_F_N_ply[1]), integrateG_phi(v_F_N_ply[2]));
		a = M * N * a_F_ply;
		break;
	case 2:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, _betaTT);
		M = 1.f;
		N = vec3(0.f);
		a = M * N * a_B_ply;
		break;
	case 3:
		M = 1.0f;
		N = vec3(0.f);
		a = M * N * AD;
		break;
	}
	return a * cosThetaI;
}

vec3 getAlphaBackward(float theta_i, float phi_i, int lobe, vec3 v_F_M_ply,
	vec3 v_F_N_ply, vec3 v_B_M_ply, vec3 a_F_ply,
	vec3 a_B_ply)
{

	float M;
	vec3 N;
	vec3 a;
	float cosThetaI = clamp(cos(theta_i), 0.f, 1.f);
	float theta_I_fiber = asin(sin(theta_i));
	switch (lobe) {
	case 0:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, fiber_bsdf->_betaR);
		M = 1.f;
		N = vec3(1.f);
		a = M * N * fR(theta_I_fiber);
		// a = a * cos(theta_I_fiber);
		break;
	case 1:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, _betaTT);
		M = 1.f;
		N = vec3(integrateG_phi(v_F_N_ply[0]), integrateG_phi(v_F_N_ply[1]), integrateG_phi(v_F_N_ply[2]));
		N = vec3(1) - N;
		a = M * N * a_F_ply;
		break;
	case 2:
		// M = integrateG(-pi / 2.0f + theta_i, pi / 2.0f + theta_i, _betaTT);
		M = 1.f;
		N = vec3(1.f);
		a = M * N * a_B_ply;
		break;
	case 3:
		M = 1.0f;
		N = vec3(1.f);
		a = M * N * AD;
		break;
	}
	return a * cosThetaI;
}

vec3 local_scattering_direct(float Num_hair_trans, Frame frame_local, vec3 wi_, vec3 wd_, float fiber_inside_, float density, vec3 trans, vec3 sigma_F, bool direct) {
	float fiber_inside = 1;
	vec3 a_F_ply;
	vec3 a_B_ply_all;
	vec3 v_F_M_ply;
	vec3 v_F_N_ply;
	float thetaI_fiber;
	vec3 v_B_M_ply_all;
	vec3 wi = to_local(frame_local, wi_);
	vec3 wd = to_local(frame_local, wd_);
	float thetaH_fiber;
	bool inside_direct;
	float thetaI;

	{
		vec3 a_R;
		vec3 a_D;
		
		vec3 wi_fiber = wi;
		{
			float sinThetaI = wi_fiber.y;
			float thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
			thetaI_fiber = thetaI;
			{
				vec3 wd_fiber = wd;
				float sinThetad = wd_fiber.y;
				float thetad_fiber = asin(clamp(sinThetad, -1.0f, 1.0f));
				thetaH_fiber = (thetad_fiber + thetaI_fiber) * 0.5f;
			}
		}
		
		{
			a_B_ply_all = vec3(0.f);
			v_B_M_ply_all = vec3(beta_diffuse);
		}
	}

	//direct bounce
	vec3 direct_bounce;
	vec3 indirect_bounce;
	{
		float SinThetaD = wd.y;
		float cosThetaD = trigInverse(SinThetaD);
		float SinThetaI = wi.y;
		float ThetaD = asin(clamp(SinThetaD, -1.f, 1.f));
		float ThetaI = asin(clamp(SinThetaI, -1.f, 1.f));
		float theta_d = (ThetaD - ThetaI) / 2;
		float theta_h = (ThetaD + ThetaI) * 0.5f;
		float phi_d = atan(wd.x, wd.z);
		float phi_i = atan(wi.x, wi.z);
		float phi = phi_d - phi_i;
		if (phi < 0.f)
			phi += 2.0f * M_PI;
		if (phi_d < 0.f)
			phi_d += 2.0f * M_PI;
		if (direct) {
			vec3 betaR = vec3(betaR);
			vec3 M_R_G = vec3(g_theta(betaR[0], thetaH_fiber), g_theta(betaR[1], thetaH_fiber),
				g_theta(betaR[2], thetaH_fiber)) *
				(fR(thetaI_fiber));
			vec3 M_D_G =
				vec3(1.0f / M_PI, 1.0f / M_PI, 1.0f / M_PI) * (AD); // forward half should be included in f Lobe
			vec3 N_R_G = vec3(0.5f / M_PI, 0.5f / M_PI, 0.5f / M_PI);
			
			vec3 N_D_G = vec3(0.5f / M_PI, 0.5f / M_PI, 0.5f / M_PI);
			vec3 beta_TT = vec3(betaTT);
			vec3 M_TT_G = vec3(g_theta(beta_TT[0], thetaH_fiber), g_theta(beta_TT[1], thetaH_fiber),
				g_theta(beta_TT[2], thetaH_fiber)) *
				(vec3(1, 1, 1) - fR(thetaI_fiber)) * ATT;
			vec3 N_TT_G = vec3(g_phi(betaN, phi - M_PI));
			
			vec3 Ret = M_R_G * N_R_G + M_TT_G * N_TT_G + M_D_G * N_D_G;
			direct_bounce = Ret;
			
		}
	}
	vec3 wd_local = to_local(frame_local, wd_);
	float sinThetaD = wd_local.y;
	float cosThetaD = trigInverse(sinThetaD);
	return (direct_bounce) * cosThetaD;
}

vec3 local_scattering_indirect(float Num_hair_trans, Frame frame_local, vec3 wi_, vec3 wd_, float fiber_inside_, float density, vec3 trans, vec3 sigma_F, bool direct) {
	float fiber_inside = 1;
	vec3 a_F_ply;
	vec3 a_B_ply_all;
	vec3 v_F_M_ply;
	vec3 v_F_N_ply;
	float thetaI_fiber;
	vec3 v_B_M_ply_all;
	vec3 wi = to_local(frame_local, wi_);
	vec3 wd = to_local(frame_local, wd_);
	float thetaH_fiber;
	bool inside_direct;
	float thetaI;
	float d_cross = 0.f;

	vec3 a_F_ply_backscattering;
	vec3 v_F_M_ply_backscattering;
	vec3 v_F_N_ply_backscattering;
	vec3 a_F_fiber;
	vec3 a_B_fiber;
	vec3 _beta_avg;
	vec3 _beta_avg_N;
	vec3 _beta_avg_backward;
	{
		vec3 a_R;
		vec3 a_D;
		vec3 b = vec3(0, 1, 0);
		vec3 ng = vec3(0, 0, 1);
		vec3 t = normalize(cross(b, ng));
		vec3 wi_fiber = wi;
		{
			float sinThetaI = wi_fiber.y;
			float thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
			thetaI_fiber = thetaI;
			{
				vec3 wd_fiber = wd;
				float sinThetad = wd_fiber.y;
				float thetad_fiber = asin(clamp(sinThetad, -1.0f, 1.0f));
				thetaH_fiber = (thetad_fiber + thetaI_fiber) * 0.5f;
			}
			vec3 AFR = getAlphaForward_Single(thetaI, 0, 0);
			vec3 AFTT = getAlphaForward_Single(thetaI, 0, 1);
			vec3 AFD = getAlphaForward_Single(thetaI, 0, 2);
			vec3 ABR = getAlphaBackward_Single(thetaI, 0, 0);
			vec3 ABTT = getAlphaBackward_Single(thetaI, 0, 1);
			vec3 ABD = getAlphaBackward_Single(thetaI, 0, 2);
			a_B_fiber = ABR + ABTT + ABD;
			a_F_fiber = AFR + AFTT + AFD;
			a_D = AFD + ABD;
			a_R = AFR + ABR;

			//  for beta_avg
			_beta_avg = (vec3(AFR) * betaR + vec3(AFTT) * betaTT +
				vec3(AFD) * beta_diffuse) /
				(vec3(AFR) + vec3(AFTT) + vec3(AFD));
			_beta_avg_backward = (vec3(ABR) * betaR + vec3(ABTT) * betaTT +
				vec3(ABD) * beta_diffuse) /
				(vec3(ABR) + vec3(ABTT) + vec3(ABD));
			_beta_avg_N =
				(vec3(AFR) * beta_diffuse + vec3(AFTT) * betaN + vec3(AFD) * beta_diffuse) /
				(vec3(AFR) + vec3(AFTT) + vec3(AFD));
			_beta_avg *= _beta_avg;
			_beta_avg_backward *= _beta_avg_backward;
			_beta_avg_N *= _beta_avg_N;       //std::cout << " thetaI_fiber is "  << thetaI_fiber << std::endl;
		}
	
	}
	//direct bounce
	
	vec3 indirect_bounce;
	//indirect bounce=
	{
		vec3 wi_local = to_local(frame_local, wi_);
		vec3 wo_local = to_local(frame_local, wd_);
		float sinThetaO = wo_local.y;
		float sinThetaI = wi_local.y;
		float theta_I = asin(clamp(sinThetaI, -1.0f, 1.0f));
		float theta_O = asin(clamp(sinThetaO, -1.0f, 1.0f));
		// Ref_bcsdf *bsdf_cur = (Ref_bcsdf *)(hit.material->bsdf);
		//  attenuation
		float theta_d = (theta_O - theta_I) / 2;
		float theta_h = (theta_O + theta_I) / 2;
		float phi_i = atan(wi_local.x, wi_local.z);
		vec3 a_B_Hat = a_B_fiber;
		vec3 a_F_Hat = a_F_fiber;
		if (length(vec3(a_F_Hat)) < 0.001) {
			indirect_bounce = vec3(0.f);
		}
		else {
			vec3 A_1 = a_B_Hat * pow(a_F_Hat, vec3(2)) / (vec3(1) - pow(a_F_Hat, vec3(2)));
			vec3 A_3 = pow(a_B_Hat, vec3(3)) * pow(a_F_Hat, vec3(2)) / (pow((vec3(1) - pow(a_F_Hat, vec3(2))), vec3(3)));
			vec3 A_B_Hat = A_1 + A_3;

			// try include TR/TTRT/TTTRTT....
			// We need to considering attenuation because of cosine term and azimuthal distribution
			vec3 A_1_flat = a_B_Hat * a_F_Hat / (vec3(1) - pow(a_F_Hat, vec3(2))) * d_cross_inter;
			vec3 A_B_Hat_flat = A_B_Hat + A_1_flat;

			// Spread
			// in our model, we don't have the shift
			float alpha_F = 0.f;
			float alpha_B = 0.f;
			float delta_B_Hat = 0.f;
			// variance
			vec3 variance_F = _beta_avg;
			vec3 variance_B = _beta_avg_backward;
			vec3 sigma_B = (vec3(1) + 0.7f * pow(a_F_Hat, vec3(2))) *
				(a_B_Hat * sqrt((2 * variance_B + variance_F)) +
					pow(a_B_Hat, vec3(3)) * sqrt((2 * variance_B + 3 * variance_F))) /
				(a_B_Hat + pow(a_B_Hat, vec3(3)) * (2 * variance_B + 3 * variance_F));
			vec3 sigma_B_flat = sqrt(sigma_B * sigma_B - vec3(variance_F));
			vec3 sigma_B_all = sigma_B * A_B_Hat / A_B_Hat_flat +
				sigma_B_flat * A_1_flat / A_B_Hat_flat;
			float phi_o = atan(wo_local.x, wo_local.z);
			float phi = phi_o - phi_i;
			if (phi < 0.f)
				phi += 2.0f * M_PI;
			if (phi_o < 0.f)
				phi_o += 2.0f * M_PI;
			if (direct) {
				vec3 S_B_Hat;
				float S_B_Hat_x = g_theta(sigma_B_all[0], theta_h);
				float S_B_Hat_y = g_theta(sigma_B_all[1], theta_h);
				float S_B_Hat_z = g_theta(sigma_B_all[2], theta_h);
				S_B_Hat = vec3(S_B_Hat_x, S_B_Hat_y, S_B_Hat_z);
				if ((phi_o > 0.5f * M_PI) && (phi_o < 1.5f * M_PI)) {
					S_B_Hat = S_B_Hat * max(0, min(cos(phi_o) + 1, 1));
				}
				vec3 local_scatter = d_b * 2.f * A_B_Hat_flat * S_B_Hat / M_PI;
				indirect_bounce = local_scatter;
			}
		}
	}
	vec3 wd_local = to_local(frame_local, wd_);
	float sinThetaD = wd_local.y;
	float cosThetaD = trigInverse(sinThetaD);
	return indirect_bounce * cosThetaD;
}

in vec4 pos;
in vec4 shadowPos;
in vec4 world_pos;
in vec3 yarn_Dir;
in vec2 uv;
uniform vec3			light_dir;
uniform vec3			view_dir;
uniform sampler2D shadow_tex;
uniform sampler2D flyaway_tex;
uniform mat4 shadow_matrix;
uniform mat4 camera_matrix;
uniform sampler2D random_tex;

uniform sampler2D env_shadow_tex_1;
uniform sampler2D env_shadow_tex_2;
uniform sampler2D env_shadow_tex_3;
uniform sampler2D env_shadow_tex_4;
uniform sampler2D env_shadow_tex_5;
uniform sampler2D env_shadow_tex_6;
uniform sampler2D env_shadow_tex_7;
uniform sampler2D env_shadow_tex_8;
uniform sampler2D env_shadow_tex_9;
uniform sampler2D env_shadow_tex_10;
uniform sampler2D env_shadow_tex_11;
uniform sampler2D env_shadow_tex_12;
uniform sampler2D env_shadow_tex_13;
uniform sampler2D env_shadow_tex_14;
uniform sampler2D env_shadow_tex_15;
uniform sampler2D env_shadow_tex_16;

uniform sampler2D integral_tex;

uniform bool use_env;
uniform vec3 envlight_dirs[16];
uniform mat4 envshadow_matrixs[16];
uniform float env_intensity[16];
uniform float env_lambdas[16];
uniform vec3 adjust_color = vec3(1,1,1);
uniform float texture_offset;
uniform vec3 residual_color[8];

layout(location = 0, index = 0) out vec4 color;

vec3 gammaCorrection (vec3 colour, float gamma) {
  return pow(colour, vec3(1. / gamma));
}

float angleBetween(vec3 v1, vec3 v2, vec3 t) {
	
    vec3 v1_proj = normalize(v1 - dot(v1, t) * t);
    vec3 v2_proj = normalize(v2 - dot(v2, t) * t);
    float dotProduct = dot(v1_proj, v2_proj);
	if(dotProduct>0.999) return 0.0;
	if(dotProduct<-0.999) return M_PI;
    float angle = acos(dotProduct);
    
    vec3 cross_product = cross(v1_proj, v2_proj);
    
    float direction = dot(cross_product, t);
    
	angle = clamp(angle, 0, M_PI);
    if(direction < 0.0) angle = M_PI_2 - angle;
	
    return angle;
}

vec3 getResidual(float angle) {
    float segment = angle * 8.0;
    int index = int(floor(segment)); // 计算当前使用的 AO 纹理对的下标
	index = index % 8;
    float mixFactor = fract(segment); // 线性插值因子
    
    vec3 color;
    
    if (index == 0) {
        color = mix(residual_color[0], residual_color[1], mixFactor);
    } else if (index == 1) {
        color = mix(residual_color[1], residual_color[2], mixFactor);
    } else if (index == 2) {
        color = mix(residual_color[2], residual_color[3], mixFactor);
    } else if (index == 3) {
        color = mix(residual_color[3], residual_color[4], mixFactor);
    } else if (index == 4) {
        color = mix(residual_color[4], residual_color[5], mixFactor);
    } else if (index == 5) {
        color = mix(residual_color[5], residual_color[6], mixFactor);
    } else if (index == 6) {
        color = mix(residual_color[6], residual_color[7], mixFactor);
    } else {
        color = mix(residual_color[7], residual_color[0], mixFactor); // circular, back to AO0
    }
    
    return color;
}

int sampleRandom(int x){
	int channel = x%4;
	int pixel = (x/4)%(1024*1024);
	int u = pixel/1024;
	int v = pixel%1024;

    vec2 uv = vec2(float(u)/1024, float(v)/1024);
    vec4 randomVec4 = texture(random_tex, uv); 
	
	float randomValue = randomVec4.w;
	if(channel==0) randomValue =randomVec4.x;
	if(channel==1) randomValue = randomVec4.y;
	if(channel==2) randomValue = randomVec4.z;

    // 将随机数 [0, 1] 映射到离散值 0, 1, 2
    int discreteValue = int(floor(randomValue * 3.0));
    discreteValue = clamp(discreteValue, 0, 2);
	return discreteValue;
}

vec3 fetch4DTexture(int invlambda, int cosphi, int costhetaI, int costhetaO){
 	int scale = 16; // 每个维度分块大小为16
    
    // 计算纹理空间的实际坐标
    int x = invlambda * scale + costhetaI;
    int y = cosphi * scale + costhetaO;

    // 提取整数部分和小数部分
    ivec2 base = ivec2(x, y);

    // 读取相邻像素的值
    vec4 c = texelFetch(integral_tex, base, 0); // 左下角值
	return c.xyz; 
}


#define lightSize 5
float PCSS(vec4 shadow_tex_coords, sampler2D shadow_tex, int kernel_size, float texel_size) {
	float shadow = 0.0;
    int half_kernel = kernel_size / 2;
	float blockerDepth = 0;
    int numBlockers = 0;
	for (int x = -half_kernel; x <= half_kernel; ++x) {
        for (int y = -half_kernel; y <= half_kernel; ++y) {
            vec4 offset_coords = shadow_tex_coords;
            offset_coords.xy += vec2(x, y) * texel_size;
			float depth_closest = texture(shadow_tex, offset_coords.xy).x;
			if((depth_closest <= offset_coords.z)){
			  numBlockers += 1;	
			  blockerDepth +=depth_closest;
		    }
        }
    }

    //if (numBlockers == 0) return .0; // 没有遮挡点，无阴影
    blockerDepth /= float(numBlockers); // 平均遮挡深度

	float currentDepth = shadow_tex_coords.z;

    float penumbraSize = (currentDepth - blockerDepth) / blockerDepth * lightSize;

	penumbraSize = 1;
    // 计算阴影模糊
    const int numPCFSamples = 25; // PCF 采样数量
    for (int x = -2; x <= 2; ++x) {
        for (int y = -2; y <= 2; ++y) {
            vec2 offset = vec2(x, y) * penumbraSize * texel_size;
            shadow += texture(shadow_tex, shadow_tex_coords.xy + offset).r >= currentDepth ? 1.0 : 0.0;
        }
    }

    return shadow / numPCFSamples;
}


float PCF(vec4 shadow_tex_coords, sampler2D shadow_tex, int kernel_size, float texel_size) {
	float shadow = 0.0;
	float currentDepth = shadow_tex_coords.z;
    const int numPCFSamples = 25; // PCF 采样数量
    for (int x = -2; x <= 2; ++x) {
        for (int y = -2; y <= 2; ++y) {
            vec2 offset = vec2(x, y) * texel_size;
            shadow += texture(shadow_tex, shadow_tex_coords.xy + offset).r >= currentDepth ? 1.0 : 0.0;
        }
    }

    return shadow / numPCFSamples;
}

vec3 sample4DTexture(float invlambda, float cosphi, float costhetaI, float costhetaO) {
    float t_invlambda = fract(invlambda);
	int invlambda_0 = clamp(int(floor(invlambda)),0,63); 
	int invlambda_1 = clamp(invlambda_0+1,0,63); 

    float t_cosphi = fract(cosphi);
	int cosphi_0 = clamp(int(floor(cosphi)),0,15); 
	int cosphi_1 = clamp(invlambda_0+1,0,15); 

    float t_costhetaI = fract(costhetaI);
	int costhetaI_0 = clamp(int(floor(costhetaI)),0,15); 
	int costhetaI_1 = clamp(costhetaI_0+1,0,15); 

    float t_costhetaO = fract(costhetaO);
	int costhetaO_0 = clamp(int(floor(costhetaO)),0,15); 
	int costhetaO_1 = clamp(costhetaO_0+1,0,15); 

	vec3 c0000 = fetch4DTexture(invlambda_0, cosphi_0, costhetaI_0, costhetaO_0);
	return c0000;
	vec3 c1000 = fetch4DTexture(invlambda_1, cosphi_0, costhetaI_0, costhetaO_0);
	vec3 c0100 = fetch4DTexture(invlambda_0, cosphi_1, costhetaI_0, costhetaO_0);
	vec3 c0010 = fetch4DTexture(invlambda_0, cosphi_0, costhetaI_1, costhetaO_0);
	vec3 c0001 = fetch4DTexture(invlambda_0, cosphi_0, costhetaI_0, costhetaO_1);
	vec3 c1100 = fetch4DTexture(invlambda_1, cosphi_1, costhetaI_0, costhetaO_0);
	vec3 c0110 = fetch4DTexture(invlambda_0, cosphi_1, costhetaI_1, costhetaO_0);
	vec3 c0011 = fetch4DTexture(invlambda_0, cosphi_0, costhetaI_1, costhetaO_1);
	vec3 c0101 = fetch4DTexture(invlambda_0, cosphi_1, costhetaI_0, costhetaO_1);
	vec3 c1010 = fetch4DTexture(invlambda_1, cosphi_0, costhetaI_1, costhetaO_0);
	vec3 c1001 = fetch4DTexture(invlambda_1, cosphi_0, costhetaI_0, costhetaO_1);
	vec3 c1110 = fetch4DTexture(invlambda_1, cosphi_1, costhetaI_1, costhetaO_0);
	vec3 c1101 = fetch4DTexture(invlambda_1, cosphi_1, costhetaI_0, costhetaO_1);
	vec3 c1011 = fetch4DTexture(invlambda_1, cosphi_0, costhetaI_1, costhetaO_1);
	vec3 c0111 = fetch4DTexture(invlambda_0, cosphi_1, costhetaI_1, costhetaO_1);
	vec3 c1111 = fetch4DTexture(invlambda_1, cosphi_1, costhetaI_1, costhetaO_1);

	vec3 cX000 = mix(c0000, c1000, t_invlambda);
	vec3 cX001 = mix(c0001, c1001, t_invlambda);
	vec3 cX010 = mix(c0010, c1010, t_invlambda);
	vec3 cX011 = mix(c0011, c1011, t_invlambda);
	vec3 cX100 = mix(c0100, c1100, t_invlambda);
	vec3 cX101 = mix(c0101, c1101, t_invlambda);
	vec3 cX110 = mix(c0110, c1110, t_invlambda);
	vec3 cX111 = mix(c0111, c1111, t_invlambda);
	 
	vec3 cXX00 = mix(cX000, cX100, t_cosphi);
	vec3 cXX01 = mix(cX001, cX101, t_cosphi);
	vec3 cXX10 = mix(cX010, cX110, t_cosphi);
	vec3 cXX11 = mix(cX011, cX111, t_cosphi);

	vec3 cXXX0 = mix(cXX00, cXX10, t_costhetaI);
	vec3 cXXX1 = mix(cXX01, cXX11, t_costhetaI);

	vec3 c = mix(cXXX0, cXXX1, t_costhetaO);
    return c;
}

float shadowCalc(vec4 world_pos, int i, int kernel_size, float texel_size){
	float global_shadow_scale;
	if(i==-1){
		vec4 world_pos_here	= shadow_matrix * world_pos;
		world_pos_here.z -= SHADOW_BIAS;
		global_shadow_scale = PCF(world_pos_here, shadow_tex, kernel_size, texel_size);
	}
	else{
		vec4 world_pos_here	= envshadow_matrixs[i] * world_pos;
		world_pos_here.z -= SHADOW_BIAS;
		if (i == 0) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_1, kernel_size, texel_size);
		if (i == 1) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_2, kernel_size, texel_size);
		if (i == 2) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_3, kernel_size, texel_size);
		if (i == 3) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_4, kernel_size, texel_size);
		if (i == 4) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_5, kernel_size, texel_size);
		if (i == 5) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_6, kernel_size, texel_size);
		if (i == 6) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_7, kernel_size, texel_size);
		if (i == 7) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_8, kernel_size, texel_size);
		if (i == 8) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_9, kernel_size, texel_size);
		if (i == 9) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_10, kernel_size, texel_size);
		if (i == 10) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_11, kernel_size, texel_size);
		if (i == 11) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_12, kernel_size, texel_size);
		if (i == 12) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_13, kernel_size, texel_size);
		if (i == 13) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_14, kernel_size, texel_size);
		if (i == 14) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_15, kernel_size, texel_size);
		if (i == 15) global_shadow_scale = PCSS(world_pos_here, env_shadow_tex_16, kernel_size, texel_size);
	}
	return global_shadow_scale;
}

float angleBetween_long(vec3 v1, vec3 t) {
	
    float dotProduct = dot(v1, t);

    float angle = acos(dotProduct);
    
	angle = clamp(angle, 0, M_PI);
	
    return angle;
}

void main()
{
	//gl_FragDepth = gl_FragCoord.z;	
	color = vec4(0,0,0,1);
	if(use_env){
		mat4 camera_matrix_inverse = inverse(camera_matrix);
		// view direction - view space
		vec4 V4 = camera_matrix * world_pos;
		vec3 V = -normalize( V4.xyz );

		// view direction - world space
		vec4 camera_pos_world = camera_matrix_inverse * vec4(0,0,0,1);
		camera_pos_world = camera_pos_world / camera_pos_world.w;
		vec3 wi1 = normalize(camera_pos_world.xyz - world_pos.xyz);
		//wi1 =vec3(0,0,1);

		///----------------------------------------------------------------------
		// diffuse - view space
		///----------------------------------------------------------------------
		mat3 cm = mat3(camera_matrix);
		mat3 cm_inverse = inverse(cm);
		vec3 T_v = normalize( cm * yarn_Dir );	
		
		int discreteValue = sampleRandom(int(floor(uv.x)));
		float fractionalPart = clamp(float(discreteValue + fract(uv.x)) / 3.0, 0.0, 1.0);
		float uvy = clamp(uv.y, 0.0, 1.0);
		vec2 real_uv = vec2(fractionalPart, uv.y);

		vec4 tex_value = texture(flyaway_tex, real_uv);
		if(tex_value.w < 0.2) discard;
		vec3 tangent_local = normalize(tex_value.xyz * 2.0f - vec3(1,1,1));

		vec3 norm_e = normalize(vec3(-T_v.y, T_v.x, 0));	//cross(z, T_v);  
		vec3 norm_V = normalize(cross(norm_e, T_v));  

		T_v = normalize(tangent_local.x * T_v + tangent_local.y * norm_e - tangent_local.z * norm_V);
		vec3 T_world = normalize(cm_inverse * T_v);

		vec3 bitangent = normalize(cross(yarn_Dir, vec3(0,0,1)));
		Frame frame_local;
		frame_local.t =  normalize(bitangent);
		frame_local.b = normalize(yarn_Dir);
		frame_local.n = normalize(cross(frame_local.t, frame_local.b));
		
		Frame frame_local_2;
		frame_local_2.t =  normalize(bitangent);
		frame_local_2.b = normalize(yarn_Dir);
		frame_local_2.n = normalize(cross(frame_local_2.t, frame_local_2.b));
		
		vec3 view_dir_local = to_local(frame_local, wi1);
		float sinThetaO = view_dir_local.y;
		float cosThetaD = trigInverse(sinThetaO);
		float thetaO = asin(clamp(sinThetaO, -1.0f, 1.0f)); 
		thetaO += 0.5 * M_PI;
		float costhetaO = cos(thetaO);
		costhetaO = (costhetaO + 1.f) * 8;
		float phi_O = atan(view_dir_local.x, view_dir_local.z);
		for(int i=0;i<4;i++)
		{
			vec3 light_dir_cur = envlight_dirs[i];
			float intensity = env_intensity[i];

			// light direction - world space
			vec3 L = normalize(light_dir_cur);

		///----------------------------- Ref_bcsdf ---------------------------------
			vec3 light_dir1 = normalize(light_dir_cur);

			vec3 light_dir_local = to_local(frame_local, light_dir1);
			vec3 fcolor;
			{
				float sinThetaI = light_dir_local.y;
			

				float thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
				thetaI += 0.5 * M_PI;
				float costhetaI = cos(thetaI);
				costhetaI = (costhetaI + 1.f) * 8;
				costhetaI = clamp(costhetaI, 0, 15);

				float phi_I = atan(light_dir_local.x, light_dir_local.z);
				float cosphi = cos(phi_I - phi_O);
				cosphi = (cosphi + 1.f)*8;
				
				float invlambda = 1 / env_lambdas[i]; 
				invlambda = invlambda * 64 - 1;
				
				fcolor = sample4DTexture(invlambda, cosphi, costhetaI, costhetaO);
				fcolor = gammaCorrection(fcolor, 1.0) * cosThetaD;
			}

			///----------------------------------------------------------------------
			//	global shadow
			///----------------------------------------------------------------------
			
			// 使用 PCF 的版本
			float texel_size = 1.0 / textureSize(shadow_tex, 0).x;  // 假设方形纹理
			int kernel_size = 7;  // 3x3 的核
			vec4 real_world_pos = shadowPos;
			
			float global_shadow_scale = shadowCalc(real_world_pos, i, kernel_size, texel_size);
			
			bool direct = true;
			vec3 trans = vec3(0,0,0);
			vec3 var = vec3(0,0,0);

			// finalize
			float fiber_frans_inside =1.f;
			float geo_unitless_density = 1.f;
			direct = true;

			if(env_lambdas[i]<3)
			fcolor = local_scattering_direct(0.f, frame_local, wi1, light_dir1, fiber_frans_inside, geo_unitless_density, trans, var, direct);
			+ local_scattering_indirect(0.f, frame_local_2, wi1, light_dir1, fiber_frans_inside, geo_unitless_density, trans, var, direct);
		
			vec3 self_coreAO_value = vec3(1,1,1);

			//global_shadow_scale = 1;
			fcolor = fcolor* self_coreAO_value * global_shadow_scale * max(0.0, intensity) * adjust_color;
			
		 	color.xyz  += fcolor;
		}
	} else
	{
		mat3 cm = mat3(camera_matrix);
		mat3 cm_inverse = inverse(cm);
		vec3 T_v = normalize( cm * yarn_Dir );			// fiber direction

		int discreteValue = sampleRandom(int(floor(uv.x)));
		float fractionalPart = clamp(float(discreteValue + fract(uv.x)) / 3.0, 0.0, 1.0);
		float uvy = clamp(uv.y, 0.0, 1.0);
		vec2 real_uv = vec2(fractionalPart, uv.y);

		vec4 tex_value = texture(flyaway_tex, real_uv);
		if(tex_value.w < 0.2) discard;
		vec3 tangent_local = normalize(tex_value.xyz * 2.0f - vec3(1,1,1));

		vec3 norm_e = normalize(vec3(-T_v.y, T_v.x, 0));	//cross(z, T_v);  
		vec3 norm_V = normalize(cross(norm_e, T_v));  

		T_v = normalize(tangent_local.x * T_v + tangent_local.y * norm_e - tangent_local.z * norm_V);
		
		vec3 T_world = normalize(cm_inverse * T_v);
	
		mat4 camera_matrix_inverse = inverse(camera_matrix);
		vec4 camera_pos_world = camera_matrix_inverse * vec4(0,0,0,1);
		camera_pos_world = camera_pos_world / camera_pos_world.w;
		vec3 wi1 = normalize(camera_pos_world.xyz - world_pos.xyz);
		//vec3 wi1 = vec3(0,0,1);
		vec3 light_dir1 = normalize(light_dir);
		//get the angle between L_v and norm_V
		vec3 L_azi = normalize(light_dir1 - dot(yarn_Dir, light_dir1) * yarn_Dir);
		vec3 face_normal = normalize(cross(normalize(cross(yarn_Dir, wi1)), yarn_Dir));
		float angle = (angleBetween(L_azi, face_normal, yarn_Dir) / (M_PI_2));
		vec3 residual = gammaCorrection(getResidual(angle), 2.18) * 2;
		
	//----------------------------- Ref_bcsdf ---------------------------------
		vec3 bitangent = normalize(cross(yarn_Dir, vec3(0,0,1)));
		Frame frame_local;
		frame_local.t =  normalize(bitangent);
		frame_local.b = normalize(yarn_Dir);
		frame_local.n = normalize(cross(frame_local.t, frame_local.b));

		Frame frame_local_2;
		frame_local_2.t =  normalize(bitangent);
		frame_local_2.b = normalize(T_world);
		frame_local_2.n = normalize(cross(frame_local_2.t, frame_local_2.b));

		float fiber_frans_inside =1.f;
		float geo_unitless_density = 1.f;
		bool direct = true;
		vec3 trans = vec3(1,1,1);
		vec3 var = vec3(0,0,0);
		vec3 fcolor = local_scattering_direct(0.f, frame_local, wi1, light_dir1, fiber_frans_inside, geo_unitless_density, trans, var, direct) * vec3(5,5,5);
		+ local_scattering_indirect(0.f, frame_local_2, wi1, light_dir1, fiber_frans_inside, geo_unitless_density, trans, var, direct) * vec3(5,5,5);
		
		float texel_size = 1.0 / textureSize(shadow_tex, 0).x;  // 假设方形纹理
		int kernel_size = 7;  // 3x3 的核
		vec4 real_world_pos = world_pos;
		
		float global_shadow_scale = shadowCalc(real_world_pos, -1, kernel_size, texel_size);

		fcolor = fcolor * global_shadow_scale * adjust_color;
		color.xyz  += fcolor;
		//color.xyz = vec3(fractionalPart, uv.y, 0);
   }
   color.xyz = gammaCorrection(color.xyz, 2.18);
   color.a = 1.0;
   //color.a = tex_value.w;
   
}
