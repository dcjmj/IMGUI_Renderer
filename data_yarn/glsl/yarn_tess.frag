#version 430 core

#define M_PI		3.141592653589793238462643383279f
#define M_PI_2		6.283185307179586476925286766559f
#define M_PI_2_INV	0.159154943091895335768883763372f

//#define HIGH_PERFORMANCE

layout(location = 0, index = 0) out vec4 color;

in GEO_OUT
{
	vec4 pos;
	vec2 uv;
	float fiber_thickness;
	vec3 offset_3d;
	vec3 ply_dir_3d;
	vec3 fbr_dir_3d;
	vec3 center_3d;	
	float reverse;
	vec4 world_pos;
	vec3 lod_color;
	vec3 up;
	float curvature;
} fs_in;

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

vec3 local_scattering(float Num_hair_trans, Frame frame_local, vec3 wi_, vec3 wd_, float fiber_inside_, float density, vec3 trans, vec3 sigma_F, bool direct) {
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
	float n_I = 0.f;
	float d_cross = 0.f;
	float density_azi = density;
	float opacity_c = 0.6f - min(density_azi, 1.0f) * 0.6f;
	float path_n = 0.f;

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
		Frame fiber_frame;
		fiber_frame.t = t;
		fiber_frame.b = b;
		fiber_frame.n = ng;
		vec3 wi_fiber = to_local(fiber_frame, wi);
		{
			float sinThetaI = wi_fiber.y;
			float thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
			thetaI_fiber = thetaI;
			{
				vec3 wd_fiber = to_local(fiber_frame, wd);
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
		
		float sinThetaI = wi.y;
		float sinThetaO = wd.y;
		float cosThetaO = trigInverse(sinThetaO);
		float cosThetaI = trigInverse(sinThetaI);
		thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
		float thetaO = asin(clamp(sinThetaO, -1.0f, 1.0f));
		float thetaD = (thetaO - thetaI) * 0.5f;
		float theta_h = (thetaI + thetaO) * 0.5f;
		// here following Marschner03's
		// float theta_h = thetaI + thetaO;
		float cosThetaD = cos(thetaD);
		float phi_I = atan(wi.x, wi.z);
		float phi_O = atan(wd.x, wd.z);
		float n_I_local;
		// l = 2h in paper, but I think it's a typo
		float phi = phi_O - phi_I;
		if (phi_O < 0.0f)
			phi_O += 2 * M_PI;
		if (phi < 0.0f)
			phi += 2 * M_PI;
		float l = 2 * cos(M_PI - phi_O);

		{
			inside_direct = true;
			n_I = 1.02f;
			n_I_local = fiber_inside;
		}
		
		// attenuation of ply
		n_I -= 1;
		n_I_local = n_I_local - 1;
		{
			a_B_ply_all = vec3(0.f);
			v_B_M_ply_all = vec3(beta_diffuse);

			a_F_ply = a_F_fiber;
			v_F_M_ply = _beta_avg;
			v_F_M_ply = min(v_F_M_ply, vec3(beta_diffuse));

			v_F_N_ply = sqrt(n_I * _beta_avg_N);
			v_F_N_ply = min(v_F_N_ply, vec3(beta_diffuse));

			float forward_fiber = max(0.0f, fiber_inside - 1);
			a_F_ply_backscattering = a_F_fiber;
			v_F_M_ply_backscattering = sqrt(_beta_avg);
			v_F_M_ply_backscattering = min(v_F_M_ply_backscattering, vec3(beta_diffuse));

			v_F_N_ply_backscattering = sqrt(_beta_avg_N);
			v_F_N_ply_backscattering = min(v_F_N_ply_backscattering, vec3(beta_diffuse));
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
			vec3 beta_B = v_B_M_ply_all;
			vec3 M_B_G =
				vec3(g_theta(beta_B[0], theta_h), g_theta(beta_B[1], theta_h), g_theta(beta_B[2], theta_h)) *
				a_B_ply_all;
			vec3 N_B_G;
			
			if (n_I < 0.1f)
				N_B_G = vec3(1.f / M_PI);
			else
				N_B_G = vec3(0.f);

			if ((phi_d >= 0.5 * M_PI) && (phi_d <= 1.5 * M_PI)) {
				float opacity;
				opacity = clamp(path_n / 1.f, 0.f, 1.f);
				vec3 Ret_backward = (opacity_c * a_F_fiber + vec3(1 - opacity_c)) *
					(M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G);
				vec3 Ret = opacity * (M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G) +
					(1.f - opacity) * Ret_backward;
				direct_bounce = Ret;
				// light_path = (opacity_c * a_F_fiber + (1 - opacity_c )) * (TT_ + R_ + D_ + B);
			}
			else {
				float opacity = max(0.f, opacity_c * float(1.f - cos(phi_d)));
				// std::cout << "phi_o is " << phi_O << std::endl;
				// std::cout << "opacity_c is " << opacity_c << std::endl;
				vec3 Ret = (opacity * a_F_fiber + vec3(1 - opacity)) *
					(M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G);
				direct_bounce = Ret;
				// light_path = 0.f;
			}
		}
		else {
			vec3 betaR = sqrt((vec3(betaR) * vec3(betaR) + sigma_F));
			vec3 M_R_G = vec3(g_theta(betaR[0], thetaH_fiber), g_theta(betaR[1], thetaH_fiber),
				g_theta(betaR[2], thetaH_fiber)) *
				(fR(thetaI_fiber));
			vec3 M_D_G =
				vec3(1.0f / M_PI, 1.0f / M_PI, 1.0f / M_PI) * (AD); // forward half should be included in f Lobe
			vec3 N_R_G = vec3(0.5f / M_PI, 0.5f / M_PI, 0.5f / M_PI);
			vec3 N_D_G = vec3(0.5f / M_PI, 0.5f / M_PI, 0.5f / M_PI);
			vec3 beta_TT = sqrt((vec3(betaTT) * vec3(betaTT) + sigma_F));
			vec3 M_TT_G = vec3(g_theta(beta_TT[0], thetaH_fiber), g_theta(beta_TT[0], thetaH_fiber),
				g_theta(beta_TT[0], thetaH_fiber)) *
				(vec3(1, 1, 1) - fR(ThetaI)) * ATT;
			vec3 N_TT_G = vec3(g_phi(betaN + sigma_F[0], phi - M_PI),
				g_phi(betaN + sigma_F[1], phi - M_PI),
				g_phi(betaN + sigma_F[2], phi - M_PI));
			vec3 beta_B = v_B_M_ply_all;
			vec3 M_B_G =
				vec3(g_theta(beta_B[0], theta_h), g_theta(beta_B[1], theta_h), g_theta(beta_B[2], theta_h)) *
				a_B_ply_all;
			vec3 N_B_G;
			if (n_I < 0.1f)
				N_B_G = vec3(1.f / M_PI);
			else
				N_B_G = vec3(0.f);
			vec3 Ret;
			if ((phi_d >= 0.5 * M_PI) && (phi_d <= 1.5 * M_PI)) {
				float opacity;
				opacity = clamp(path_n / 1.f, 0.f, 1.f);
				vec3 Ret_backward = (opacity_c * a_F_fiber + vec3(1 - opacity_c)) *
					(M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G);
				Ret = opacity * (M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G) +
					(1.f - opacity) * Ret_backward;
				// light_path = (opacity_c * a_F_fiber + (1 - opacity_c )) * (TT_ + R_ + D_ + B);
			}
			else {
				float opacity = max(0.f, opacity_c * float(1.f - cos(phi_d)));
				// std::cout << "phi_o is " << phi_O << std::endl;
				// std::cout << "opacity_c is " << opacity_c << std::endl;
				Ret = (opacity * a_F_fiber + vec3(1 - opacity)) *
					(M_R_G * N_R_G + M_TT_G * N_TT_G + M_B_G * N_B_G + M_D_G * N_D_G);
				// light_path = 0.f;
			}
			Ret = trans * d_f * Ret;
			direct_bounce = Ret;
		}
	}
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
			else {
				sigma_B_all = sqrt((sigma_B_all * sigma_B_all + sigma_F));
				vec3 S_B_Hat;
				float S_B_Hat_x = g_theta(sigma_B_all[0], theta_h);
				float S_B_Hat_y = g_theta(sigma_B_all[1], theta_h);
				float S_B_Hat_z = g_theta(sigma_B_all[2], theta_h);
				S_B_Hat = vec3(S_B_Hat_x, S_B_Hat_y, S_B_Hat_z);
				if ((phi_o > 0.5f * M_PI) && (phi_o < 1.5f * M_PI)) {
					S_B_Hat = S_B_Hat * max(0, min(cos(phi_o) + 1, 1));
				}
				vec3 local_scatter = d_b * 2.f * A_B_Hat_flat * S_B_Hat / M_PI;
				// do I need the M_PI?
				// local_scatter = 1.0f * M_PI * local_scatter;
				//ASSERT(local_scatter.allFinite() && (local_scatter >= 0.0f).all());
				indirect_bounce = d_f * trans * local_scatter * M_PI;
			}
		}
	}
	vec3 wd_local = to_local(frame_local, wd_);
	float sinThetaD = wd_local.y;
	float cosThetaD = trigInverse(sinThetaD);
	//return direct_bounce* cosThetaD;
	//return indirect_bounce* cosThetaD;
	return (direct_bounce + indirect_bounce) * cosThetaD;
}

uniform bool			use_ao;
uniform bool			use_diffuse;
uniform bool			use_specular;
uniform bool			use_shadow;
uniform bool			use_self_shadow;
uniform bool			use_lod_vis;
uniform float			scale;

uniform vec3			light_dir;
uniform vec3			light_pos;
uniform vec3			light_pos_world;
uniform vec3			view_dir;

uniform mat4			shadow_matrix;
uniform mat4			view_matrix;
uniform mat4			camera_matrix;
uniform mat4			inv_view_matrix;
uniform mat4			shadow_R_matrix;

uniform sampler2D		shadow_R_tex;
uniform sampler2D		core_tex;
uniform sampler2D		fbr_dir_tex;
uniform sampler2DShadow shadow_tex;

uniform sampler3D		Core_AO_tex;
//---------------------------------       azimuthal Core AO ---------------------------------
uniform sampler2D		Core_AO0_tex;
uniform sampler2D		Core_AO1_tex;
uniform sampler2D		Core_AO2_tex;
uniform sampler2D		Core_AO3_tex;
uniform sampler2D		Core_AO4_tex;
uniform sampler2D		Core_AO5_tex;
uniform sampler2D		Core_AO6_tex;
uniform sampler2D		Core_AO7_tex;


//---------------------------------       longitudinal Core AO ---------------------------------
uniform sampler2D		Core_AO0_tex_Long;
uniform sampler2D		Core_AO1_tex_Long;
uniform sampler2D		Core_AO2_tex_Long;
uniform sampler2D		Core_AO3_tex_Long;
uniform sampler2D		Core_AO4_tex_Long;
uniform sampler2D		Core_AO5_tex_Long;
uniform sampler2D		Core_AO6_tex_Long;
uniform sampler2D		Core_AO7_tex_Long;
uniform sampler2D		Core_AO8_tex_Long;

uniform vec3			default_color	= vec3(0.93, 0.53, 0.01);
uniform float			lamda_R			= 0.01;
uniform float			lamda_TT		= 0.28;
uniform float			k_R				= 0.1;
uniform float			k_TT			= 0.001;
uniform float mDelta_0;
uniform float mDelta_1;
uniform float mDelta_2;

vec3 V, L;
vec3 T_v, L_v, N_v;
//vec3 z = vec3(0,0,1);
vec3 defaultColor = vec3(0.2f, 0.8f,0.1f);
float max_fiber_radius = 0.07f;
float max_fiber_radius_inv = 14.28f;

float diffuse, specular;
//#define KAJIYA_KAY


vec3 gammaCorrection (vec3 colour, float gamma) {
  return pow(colour, vec3(1. / gamma));
}

float angleBetween(vec3 v1, vec3 v2, vec3 t) {
	
    vec3 v1_proj = v1 - dot(v1, t) * t;
    vec3 v2_proj = v2 - dot(v2, t) * t;
    float dotProduct = dot(v1_proj, v2_proj);

    float angle = acos(dotProduct);
    
    vec3 cross_product = cross(v1_proj, v2_proj);
    
    float direction = dot(cross_product, t);
    
	angle = clamp(angle, 0, M_PI);
    if(direction < 0.0) angle = M_PI_2 - angle;
	
    return angle;
}

float angleBetween_long(vec3 v1, vec3 t) {
	
    float dotProduct = dot(v1, t);

    float angle = acos(dotProduct);
    
	angle = clamp(angle, 0, M_PI);
	
    return angle;
}

vec4 getAOTexture(sampler2D AO0, sampler2D AO1, sampler2D AO2, sampler2D AO3,
                  sampler2D AO4, sampler2D AO5, sampler2D AO6, sampler2D AO7,
                  vec2 uv, float angle) {
    float segment = angle * 8.0;
    int index = int(floor(segment)); // 计算当前使用的 AO 纹理对的下标
    float mixFactor = fract(segment); // 线性插值因子
    
    vec4 color;
    
    if (index == 0) {
        color = mix(texture(AO0, uv), texture(AO1, uv), mixFactor);
    } else if (index == 1) {
        color = mix(texture(AO1, uv), texture(AO2, uv), mixFactor);
    } else if (index == 2) {
        color = mix(texture(AO2, uv), texture(AO3, uv), mixFactor);
    } else if (index == 3) {
        color = mix(texture(AO3, uv), texture(AO4, uv), mixFactor);
    } else if (index == 4) {
        color = mix(texture(AO4, uv), texture(AO5, uv), mixFactor);
    } else if (index == 5) {
        color = mix(texture(AO5, uv), texture(AO6, uv), mixFactor);
    } else if (index == 6) {
        color = mix(texture(AO6, uv), texture(AO7, uv), mixFactor);
    } else {
        color = mix(texture(AO7, uv), texture(AO0, uv), mixFactor); // circular, back to AO0
    }
    
    return color;
}

vec4 getAOTexture_Long(sampler2D AO0, sampler2D AO1, sampler2D AO2, sampler2D AO3,
                  sampler2D AO4, sampler2D AO5, sampler2D AO6, sampler2D AO7, sampler2D AO8,
                  vec2 uv, float angle) {
    float segment = angle * 8.0;
    int index = int(floor(segment)); // 计算当前使用的 AO 纹理对的下标
    float mixFactor = fract(segment); // 线性插值因子
    
    vec4 color;
    
    if (index == 0) {
        color = mix(texture(AO0, uv), texture(AO1, uv), mixFactor);
    } else if (index == 1) {
        color = mix(texture(AO1, uv), texture(AO2, uv), mixFactor);
    } else if (index == 2) {
        color = mix(texture(AO2, uv), texture(AO3, uv), mixFactor);
    } else if (index == 3) {
        color = mix(texture(AO3, uv), texture(AO4, uv), mixFactor);
    } else if (index == 4) {
        color = mix(texture(AO4, uv), texture(AO5, uv), mixFactor);
    } else if (index == 5) {
        color = mix(texture(AO5, uv), texture(AO6, uv), mixFactor);
    } else if (index == 6) {
        color = mix(texture(AO6, uv), texture(AO7, uv), mixFactor);
    } else {
        color = mix(texture(AO7, uv), texture(AO8, uv), mixFactor); 
    }
    
    return color;
}

float PCF(vec4 shadow_tex_coords, sampler2DShadow shadow_tex, int kernel_size, float texel_size) {
    float shadow = 0.0;
    int half_kernel = kernel_size / 2;
    for (int x = -half_kernel; x <= half_kernel; ++x) {
        for (int y = -half_kernel; y <= half_kernel; ++y) {
            vec4 offset_coords = shadow_tex_coords;
            offset_coords.xy += vec2(x, y) * shadow_tex_coords.w * texel_size;
            shadow += textureProj(shadow_tex, offset_coords).r;
        }
    }
    return shadow / (kernel_size * kernel_size);
}

float curve_up_exponential(float x) {
    float factor = 2.0;  // 控制曲线的陡峭程度
    return 1.0 - exp(-factor * x);
}

#define SHADOW_BIAS 0.0002f      

void main()
{
	mat4 camera_matrix_inverse = inverse(camera_matrix);


	// fetch texture
	vec4 tex_value = texture(core_tex, fs_in.uv);
	vec4 tangent_tex_value = texture(fbr_dir_tex, fs_in.uv);
	vec3 tangent_local = normalize( tangent_tex_value.xyz * 2.0f - vec3(1));
	float alpha_T = tangent_tex_value.w;

	if (tex_value.w < 0.5)	discard;

	// light direction - world space
	L = normalize(light_dir);

	// view direction - view space
	vec4 V4 = camera_matrix * vec4(fs_in.center_3d,1);
	V = -normalize( V4.xyz );

	// view direction - world space
	vec4 camera_pos_world = camera_matrix_inverse * vec4(0,0,0,1);
	camera_pos_world = camera_pos_world / camera_pos_world.w;
	vec3 wi1 = normalize(camera_pos_world.xyz - fs_in.center_3d);
	wi1 =vec3(0,0,1);

	///----------------------------------------------------------------------
	// diffuse - view space
	///----------------------------------------------------------------------
	mat3 cm = mat3(camera_matrix);
	mat3 cm_inverse = inverse(cm);
	T_v = normalize( cm * fs_in.fbr_dir_3d );			// fiber direction
	L_v = normalize( cm * L );

	N_v = normalize( tex_value.xyz * 2.0f - vec3(1));
	vec3 norm_e = normalize(vec3(-T_v.y, T_v.x, 0));	//cross(z, T_v);  
	vec3 norm_V = normalize(cross(norm_e, T_v));  

	//get the angle between L_v and norm_V
	vec3 L_azi = normalize(L - dot(fs_in.fbr_dir_3d, L) * fs_in.fbr_dir_3d);
	vec3 face_normal = normalize(cross(normalize(cross(fs_in.fbr_dir_3d, wi1)), fs_in.fbr_dir_3d));
	float angle = (angleBetween(L_azi, face_normal, fs_in.fbr_dir_3d) / (M_PI_2));

	//vec3 L_long = normalize(L_v - dot(norm_V, L_v) * norm_V);
	float angle_long = 1 - (angleBetween_long(L, fs_in.fbr_dir_3d) / (M_PI));

	vec3 self_coreAO_Azi_value = getAOTexture(Core_AO0_tex, Core_AO1_tex, Core_AO2_tex, Core_AO3_tex,
                  Core_AO4_tex, Core_AO5_tex, Core_AO6_tex, Core_AO7_tex,
                  fs_in.uv, angle).xyz;
	self_coreAO_Azi_value = gammaCorrection(self_coreAO_Azi_value, 1.0) * 2;

	vec3 self_coreAO_Long_value = getAOTexture_Long(Core_AO0_tex_Long, Core_AO1_tex_Long, Core_AO2_tex_Long, Core_AO3_tex_Long,
                  Core_AO4_tex_Long, Core_AO5_tex_Long, Core_AO6_tex_Long, Core_AO7_tex_Long, Core_AO8_tex_Long,
                  fs_in.uv, angle_long).xyz;
	self_coreAO_Long_value =  gammaCorrection(self_coreAO_Long_value, 1.0) * 2;

	N_v = normalize(N_v.x *T_v + N_v.y * norm_e - N_v.z * norm_V);
	vec4 N_world_v4 = camera_matrix_inverse * vec4(N_v, 1.0);
	vec3 N_world = N_world_v4.xyz / N_world_v4.w;
	N_world = normalize(cm_inverse * N_v);

	T_v = normalize(tangent_local.x * T_v + tangent_local.y * norm_e - tangent_local.z * norm_V);
	vec3 T_world = normalize(cm_inverse * T_v);
	diffuse = max(dot(N_v, L_v), 0.0);
	
	///----------------------------------------------------------------------
	// specular
	///----------------------------------------------------------------------
	vec3 L_prime = normalize( L_v - dot(L_v, T_v) );
	vec3 V_prime = normalize( V - dot(V, T_v) );
	/*
	float a = acos(dot(-L_v, T_v)) - acos(dot(V, T_v));
	float M = exp(- a * a * 0.5 / lamda_R / lamda_R) / sqrt(M_PI_2) / lamda_R;
	specular = max(k_R * M * dot(L_prime, V_prime), 0.0);
	*/

	float phi = acos(dot(L_prime, V_prime));
	float theta_h = ( acos(dot(L_v, T_v)) + acos(dot(V, T_v)) ) * 0.5 - 1.57;
	specular = max(k_R * cos(phi * 0.5) * exp(- theta_h * theta_h * 0.5 / (lamda_R * lamda_R) ), 0.0);


#ifdef KAJIYA_KAY
	float c = dot(T_v, L_v);
	float s = sqrt(1 - c * c);
	vec3 H = normalize(V + L_v);
	float sc = dot(T_v, H);
	float ss = sqrt(1 - sc * sc);
	diffuse = s;
	specular = s * pow(ss, 38.0) * 0.6;
#endif

///----------------------------- Ref_bcsdf ---------------------------------
	//T_world = fs_in.fbr_dir_3d;
	vec3 light_dir1 = normalize(light_dir);
	vec3 bitangent = normalize(cross(T_world, N_world));
	Frame frame_local;
	frame_local.t = bitangent;
	frame_local.b = normalize(T_world);
	frame_local.n = normalize(N_world);

///----------------------------------------------------------------------------

	///----------------------------------------------------------------------
	//	global shadow
	///----------------------------------------------------------------------
	
	// 使用 PCF 的版本
	float texel_size = 1.0 / textureSize(shadow_tex, 0).x;  // 假设方形纹理
	int kernel_size = 7;  // 3x3 的核
	vec4 real_world_pos = fs_in.world_pos;
	real_world_pos.xyz += fs_in.fiber_thickness * abs(0.5 - fs_in.uv.y) * light_dir1;
	vec4 world_pos_here;
	world_pos_here	= shadow_matrix * real_world_pos;

	world_pos_here.z -= SHADOW_BIAS;

	float global_shadow_scale = PCF(world_pos_here, shadow_tex, kernel_size, texel_size);


	///----------------------------------------------------------------------
	//	global scattering
	///----------------------------------------------------------------------
	//forward_trans
	vec3 one_trans;
	vec3 one_var;
	{
		vec3 wd_ = light_dir1;
		vec3 wd = to_local(frame_local, wd_);
		vec3 wi = to_local(frame_local, wi1);
		vec3 a_B_fiber;
		vec3 a_F_fiber;
		vec3 a_R;
		vec3 a_D;
		vec3 _beta_avg;
		vec3 _beta_avg_N;
		vec3 _beta_avg_backward;
		vec3 b = vec3(0, 1, 0);
		vec3 ng = vec3(0, 0, 1);
		vec3 t = normalize(cross(b, ng));
		Frame fiber_frame;
		fiber_frame.t = t;
		fiber_frame.b = b;
		fiber_frame.n = ng;
		vec3 wd_fiber = to_local(frame_local, wd);
		{
			float sinThetaI = wd_fiber.y;
			float thetaI = asin(clamp(sinThetaI, -1.0f, 1.0f));
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

			one_trans = (AFR + AFTT + AFD);

			// variance of forward scattering = sum of avg variance of each event

			one_var = (AFR * betaR + AFTT * betaTT +
				AFD * beta_diffuse) /
				one_trans;
			one_var *= one_var;
		}
	}		

	bool direct = true;
	vec3 trans = vec3(0,0,0);
	vec3 var = vec3(0,0,0);

	///----------------------------------------------------------------------
	//	KAJIYA_KAY
	///----------------------------------------------------------------------
	vec3 color_c = default_color;

	///----------------------------------------------------------------------
	// DISABLE HERE TO BE FASTER !!!
	//*
#ifndef HIGH_PERFORMANCE
	if(!use_diffuse) diffuse = 1.0f;
	if(!use_shadow) global_shadow_scale = 1.0f;
	if(!use_specular) specular = 0.0;
	if(use_lod_vis) color_c = fs_in.lod_color;
#endif
	//*/

	// finalize
	float fiber_frans_inside =1.f;
	float geo_unitless_density = 1.f;
	direct = true;
	vec3 fcolor = local_scattering(0.f, frame_local, wi1, light_dir1, fiber_frans_inside, geo_unitless_density, trans, var, direct) * vec3(8,8,8);

	vec3 self_coreAO_value = self_coreAO_Long_value * self_coreAO_Azi_value;
	//self_coreAO_value = self_coreAO_Azi_value;

	color.xyz = fcolor * self_coreAO_value * global_shadow_scale * 0.8 + AD  * 0.2;
	color.xyz = fcolor * self_coreAO_value * global_shadow_scale * 0.8 + AD  * self_coreAO_value * 0.2;
	color.xyz = fcolor * self_coreAO_value * global_shadow_scale * 0.7 + AD  * self_coreAO_value * 0.3;
	color.xyz = vec3(self_coreAO_value);
	color.xyz = fcolor * self_coreAO_value * global_shadow_scale * mDelta_1 + AD  * self_coreAO_value * (1-mDelta_1);
	vec3 M_D_G =
		vec3(1.0f / M_PI, 1.0f / M_PI, 1.0f / M_PI) * (AD); // forward half should be included in f Lobe
	
	vec3 N_D_G = vec3(0.5f / M_PI, 0.5f / M_PI, 0.5f / M_PI);
	vec3 D = M_D_G * N_D_G;
	color.xyz = fcolor * self_coreAO_value * global_shadow_scale * mDelta_1;
	//float angle = 1 - (angleBetween(L_azi, face_normal, fs_in.fbr_dir_3d) / (M_PI_2));	
	color.xyz = gammaCorrection(color.xyz, 2.18);

	//color.xyz = fs_in.lod_color;
	color.a = 1.0;
}
