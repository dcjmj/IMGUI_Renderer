#pragma once

#include <omp.h>
#include <Eigen/Geometry>
#include "Rng.h"
#include <iostream>
#include <cstdio>

#define DEG2RAD(v) ((v)*float(.017453292519943295f))
extern std::string WORK_PATH;
namespace ks
{

	template <int N>
	using vec = Eigen::Matrix<float, N, 1>;
	using vec2 = Eigen::Vector2f;
	using vec2i = Eigen::Vector2i;
	using vec3 = Eigen::Vector3f;
	using vec3i = Eigen::Vector3i;
	using vec4 = Eigen::Vector4f;
	using vec4i = Eigen::Vector4i;
	using vec3d = Eigen::Vector3d;
	using mat2 = Eigen::Matrix2f;
	using mat3 = Eigen::Matrix3f;
	using mat4 = Eigen::Matrix4f;
	using quat = Eigen::Quaternionf;
	using quatd = Eigen::Quaterniond;
	template <int N>
	using color = Eigen::Array<float, N, 1>;
	using color3 = color<3>;
	using color4 = color<4>;
	template <int N>
	using colord = Eigen::Array<double, N, 1>;
	using color3d = colord<3>;
	using color4d = colord<4>;
	template <int N>
	using arr = Eigen::Array<float, N, 1>;
	using arr2 = arr<2>;
	using arr3 = arr<3>;
	using arr4 = arr<4>;
	template <int N>
	using arri = Eigen::Array<int, N, 1>;
	using arr2i = arri<2>;
	using arr3i = arri<3>;
	using arr4i = arri<4>;
	template <int N>
	using arru = Eigen::Array<uint32_t, N, 1>;
	using arr2u = arru<2>;
	using arr3u = arru<3>;
	using arr4u = arru<4>;

	inline constexpr float e = 2.71828182845904523536f;
	inline constexpr float pi = 3.14159265358979323846f;
	inline constexpr float beta_diffuse = pi / 3.464f;
	inline constexpr float two_pi = 6.28318530718f;
	inline constexpr float half_pi = pi / 2.0f;
	inline constexpr float quarter_pi = pi / 4.0f;
	inline constexpr float sqrt_2 = 1.41421356237f;
	inline constexpr float inv_pi = 1.0f / pi;
	inline constexpr float inf = std::numeric_limits<float>::infinity();
	inline const float before_one = std::nextafter(1.0f, -std::numeric_limits<float>::infinity());
	inline constexpr int int_max = std::numeric_limits<int>::max();

	const float eps = 1e-5f;
	const unsigned long long giga = 1024ULL * 1024ULL * 1024ULL;

	static std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems) {
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
		return elems;
	}


	static std::vector<std::string> split(const std::string& s, char delim) {
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	}

	//extern std::vector<Rng> rngs;

	static void wait(int seconds) {
		int endwait;
		endwait = clock() + seconds * CLOCKS_PER_SEC;
		while (clock() < endwait) {
			// NO-OP
		}
	}

	static float rand01() {
		int threadID = omp_get_thread_num();
		//assert(threadID >= 0 && threadID < rngs.size());
#if 0
		return rngs[threadID].rand(0.f, 1.f);
#else
	//srand(time(NULL));
		return (float)((double)rand() / (RAND_MAX));
#endif
	}

	template <typename T>
	constexpr T sqr(const T& x)
	{
		return x * x;
	}

	// NOTE: this returns {-1, 0, 1}.
	template <typename T>
	constexpr T sgn(const T& x)
	{
		return T((T(0) < x) - (x < T(0)));
	}

	// NOTE: this only returns {-1, 1}. For +/-0 it also returns the FP sign.
	template <typename T>
	T signum(const T& x)
	{
		return std::copysign(T(1.0), x);
	}

	template <typename T>
	constexpr T saturate(const T& x)
	{
		return std::clamp(x, T(0), T(1));
	}

	template <typename T>
	T safe_sqrt(const T& x)
	{
		return std::sqrt(std::max(T(0), x));
	}

	template <typename T>
	T safe_acos(T x)
	{
		return std::acos(std::clamp(x, T(-1.0), T(1.0)));
	}

	template <typename T>
	inline T mod(T a, T b)
	{
		T result = a - std::floor(a / b) * b;
		return (T)((result < 0) ? result + b : result);
	}

	// floor division: divide and round toward -inf
	// https://stackoverflow.com/questions/3041946/how-to-integer-divide-round-negative-numbers-down
	// works for any b != 0
	constexpr int divF(int a, int b) { return a / b - (sgn(a % b) == -sgn(b)); }

	template <typename T>
	inline T fract(T x)
	{
		return x - std::floor(x);
	}

	// template <typename DerivedV, typename DerivedB>
	// auto clamp(const Eigen::ArrayBase<DerivedV>& v, const Eigen::ArrayBase<DerivedB>& low,
	//     const Eigen::ArrayBase<DerivedB>& high)
	// {
	//     return v.min(high).max(low);
	// }

	template <typename T>
	inline T clamp(const T& val, const T& low, const T& high)
	{
		return std::clamp(val, low, high);
	}

	template <typename Derived>
	inline auto clamp_negative(const Eigen::ArrayBase<Derived>& v)
	{
		return v.max(v.Zero());
	}

	template <typename T>
	inline T clamp_negative(const T& v)
	{
		return std::max(v, T(0));
	}

	// Is this stable?
	template <typename T>
	inline T sinc(T x)
	{
		if (x == T(0)) {
			return T(1);
		}
		using std::sin;
		return sin(x) / x;
	}

	template <typename T>
	constexpr bool is_pow2(T v)
	{
		return v && !(v & (v - 1));
	}

	constexpr int32_t round_up_pow2(int32_t v)
	{
		v--;
		v |= v >> 1;
		v |= v >> 2;
		v |= v >> 4;
		v |= v >> 8;
		v |= v >> 16;
		return v + 1;
	}

	inline int log2_int(uint32_t v)
	{
#if defined(_MSC_VER)
		unsigned long lz = 0;
		if (_BitScanReverse(&lz, v))
			return lz;
		return 0;
#else
		return 31 - __builtin_clz(v);
#endif
	}

	inline int log2_int(int32_t v) { return log2_int((uint32_t)v); }

	inline float int_as_float(int i)
	{
		float f;
		std::memcpy(&f, &i, 4);
		return f;
	}

	inline int float_as_int(float f)
	{
		int i;
		std::memcpy(&i, &f, 4);
		return i;
	}

	inline float uint_as_float(uint32_t u)
	{
		float f;
		std::memcpy(&f, &u, 4);
		return f;
	}

	inline uint32_t float_as_uint(float f)
	{
		uint32_t u;
		std::memcpy(&u, &f, 4);
		return u;
	}

	// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
	// Also pbrt hair doc
	// TODO: as of Haswell, the PEXT instruction could do all this in a
	// single instruction.

	// "Insert" a 0 bit after each of the 16 low bits of x
	constexpr uint32_t part_1_by_1(uint32_t x)
	{
		x &= 0x0000ffff;                 // x = ---- ---- ---- ---- fedc ba98 7654 3210
		x = (x ^ (x << 8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
		x = (x ^ (x << 4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
		x = (x ^ (x << 2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
		x = (x ^ (x << 1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
		return x;
	}

	// "Insert" two 0 bits after each of the 10 low bits of x
	constexpr uint32_t part_1_by_2(uint32_t x)
	{
		x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
		x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x << 8)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x << 4)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x << 2)) & 0x09249249;  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		return x;
	}

	// Inverse of part_1_by_1 - "delete" all odd-indexed bits
	constexpr uint32_t compact_1_by_1(uint32_t x)
	{
		x &= 0x55555555;                 // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
		x = (x ^ (x >> 1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
		x = (x ^ (x >> 2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
		x = (x ^ (x >> 4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
		x = (x ^ (x >> 8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
		return x;
	}

	// Inverse of part_1_by_2 - "delete" all bits not at positions divisible by 3
	constexpr uint32_t compact_1_by_2(uint32_t x)
	{
		x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		x = (x ^ (x >> 2)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x >> 4)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x >> 8)) & 0xff0000ff;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
		return x;
	}

	constexpr uint32_t encode_morton_2(uint32_t x, uint32_t y) { return (part_1_by_1(y) << 1) + part_1_by_1(x); }

	constexpr uint32_t encode_morton_3(uint32_t x, uint32_t y, uint32_t z)
	{
		return (part_1_by_2(z) << 2) + (part_1_by_2(y) << 1) + part_1_by_2(x);
	}

	constexpr void decode_morton_2(uint32_t code, uint32_t& x, uint32_t& y)
	{
		x = compact_1_by_1(code >> 0);
		y = compact_1_by_1(code >> 1);
	}

	constexpr void decode_morton_3(uint32_t code, uint32_t& x, uint32_t& y, uint32_t& z)
	{
		x = compact_1_by_2(code >> 0);
		y = compact_1_by_2(code >> 1);
		z = compact_1_by_2(code >> 2);
	}

	inline vec2 demux_float(float f)
	{
		// ASSERT(f >= 0 && f < 1);
		uint64_t v = f * (1ull << 32);
		// ASSERT(v < 0x100000000);
		uint32_t bits[2] = { compact_1_by_1(v), compact_1_by_1(v >> 1) };
		return { bits[0] / float(1 << 16), bits[1] / float(1 << 16) };
	}

	constexpr float to_radian(float degree) { return degree / 180.0f * pi; }

	constexpr float to_degree(float radian) { return radian / pi * 180.0f; }

	inline void to_spherical(const vec3& dir, float& phi, float& theta)
	{
		theta = std::acos(std::clamp(dir.z(), -1.0f, 1.0f));
		phi = std::atan2(dir.y(), dir.x());
		if (phi < 0.0f) {
			phi += two_pi;
		}
	}

	inline vec3 to_cartesian(float phi, float theta)
	{
		float cos_theta = std::cos(theta);
		float sin_theta = std::sin(theta);
		float cos_phi = std::cos(phi);
		float sin_phi = std::sin(phi);
		return vec3(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta);
	}

	// Sometimes I want the up axis to be +y
	inline void to_spherical_yup(const vec3& dir, float& phi, float& theta)
	{
		theta = std::acos(std::clamp(dir.y(), -1.0f, 1.0f));
		phi = std::atan2(dir.z(), dir.x());
		if (phi < 0.0f) {
			phi += two_pi;
		}
	}

	inline vec3 to_cartesian_yup(float phi, float theta)
	{
		float cos_theta = std::cos(theta);
		float sin_theta = std::sin(theta);
		float cos_phi = std::cos(phi);
		float sin_phi = std::sin(phi);
		return vec3(cos_phi * sin_theta, cos_theta, sin_phi * sin_theta);
	}

	inline vec3 square_to_hemisphere(vec2 u)
	{
		// Map uniform random numbers to $[-1,1]^2$
		u = 2.0f * u - vec2::Constant(1.0f);

		// Handle degeneracy at the origin
		if (u.x() == 0.0f && u.y() == 0.0f) {
			return vec3(0.0f, 0.0f, 1.0f);
		}

		// Apply concentric mapping to point
		float phi, r;
		if (std::abs(u.x()) > std::abs(u.y())) {
			r = u.x();
			phi = (pi * 0.25f) * (u.y() / u.x());
		}
		else {
			r = u.y();
			phi = (pi * 0.5f) - (pi * 0.25f) * (u.x() / u.y());
		}

		float r2 = sqr(r);
		float sin_theta = r * std::sqrt(2.0f - sqr(r));

		float x = std::cos(phi) * sin_theta;
		float y = std::sin(phi) * sin_theta;
		float z = 1.0f - r2;
		return vec3(x, y, z);
	}

	inline vec2 hemisphere_to_square(const vec3& w)
	{
		float r = safe_sqrt(1.0f - w.z());
		float phi = std::atan2(w.y(), w.x());
		if (phi < -0.25f * pi) {
			phi += two_pi;
		}
		float x, y;
		if (phi < 0.25f * pi) {
			x = r;
			y = 4.0f / pi * r * phi;
		}
		else if (phi < 0.75f * pi) {
			x = -4.0f / pi * r * (phi - 0.5f * pi);
			y = r;
		}
		else if (phi < 1.25f * pi) {
			x = -r;
			y = -4.0f / pi * r * (phi - pi);
		}
		else {
			x = 4.0f / pi * r * (phi - 1.5f * pi);
			y = -r;
		}
		x = x * 0.5f + 0.5f;
		y = y * 0.5f + 0.5f;
		return vec2(x, y);
	}

	inline float delta_angle(float alpha, float beta)
	{
		float phi = std::fmod(std::abs(beta - alpha), two_pi);
		float distance = phi > pi ? two_pi - phi : phi;
		return distance;
	}

	inline bool solve_linear_system_2x2(const float A[2][2], const float B[2], float& x0, float& x1)
	{
		float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
		if (std::abs(det) < 1e-10f)
			return false;
		x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
		x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
		if (std::isnan(x0) || std::isnan(x1))
			return false;
		return true;
	}

	inline bool solve_quadratic(float a, float b, float c, float& t0, float& t1)
	{
		if (a == 0.0f && b == 0.0f) {
			return false;
		}
		if (a == 0.0f) {
			t0 = t1 = -c / b;
		}
		// Find quadratic discriminant
		double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
		if (discrim < 0)
			return false;
		discrim = std::sqrt(discrim);

		// Compute quadratic _t_ values
		double q;
		if (b < 0)
			q = -.5 * (b - discrim);
		else
			q = -.5 * (b + discrim);
		t0 = q / a;
		t1 = c / q;
		if (t0 > t1)
			std::swap(t0, t1);
		return true;
	}

	inline float erfinv(float x)
	{
		float w, p;
		x = std::clamp(x, -.99999f, .99999f);
		w = -std::log((1 - x) * (1 + x));
		if (w < 5) {
			w = w - 2.5f;
			p = 2.81022636e-08f;
			p = 3.43273939e-07f + p * w;
			p = -3.5233877e-06f + p * w;
			p = -4.39150654e-06f + p * w;
			p = 0.00021858087f + p * w;
			p = -0.00125372503f + p * w;
			p = -0.00417768164f + p * w;
			p = 0.246640727f + p * w;
			p = 1.50140941f + p * w;
		}
		else {
			w = std::sqrt(w) - 3;
			p = -0.000200214257f;
			p = 0.000100950558f + p * w;
			p = 0.00134934322f + p * w;
			p = -0.00367342844f + p * w;
			p = 0.00573950773f + p * w;
			p = -0.0076224613f + p * w;
			p = 0.00943887047f + p * w;
			p = 1.0167406f + p * w;
			p = 2.83297682f + p * w;
		}
		return p * x;
	}

	inline float cross2(const vec2& a, const vec2& b) { return a.x() * b.y() - a.y() * b.x(); }

	constexpr float srgb_to_linear(float x)
	{
		if (x < 0.04045f) {
			return x / 12.92f;
		}
		else {
			return std::pow((x + 0.055f) / 1.055f, 2.4f);
		}
	}

	inline float luminance(const color3& rgb)
	{
		constexpr float lum_weight[3] = { 0.212671f, 0.715160f, 0.072169f };
		return lum_weight[0] * rgb[0] + lum_weight[1] * rgb[1] + lum_weight[2] * rgb[2];
	}

	template <typename DerivedA, typename DerivedB>
	auto lerp(const Eigen::MatrixBase<DerivedA>& v1, const Eigen::MatrixBase<DerivedB>& v2, float t)
	{
		return (1.0f - t) * v1 + t * v2;
	}

	template <typename DerivedA, typename DerivedB>
	auto lerp(const Eigen::ArrayBase<DerivedA>& v1, const Eigen::ArrayBase<DerivedB>& v2, float t)
	{
		return (1.0f - t) * v1 + t * v2;
	}
	// We can have per-channel lerp for arrays
	template <typename DerivedA, typename DerivedB, typename DerivedC>
	auto lerp(const Eigen::ArrayBase<DerivedA>& v1, const Eigen::ArrayBase<DerivedB>& v2,
		const Eigen::ArrayBase<DerivedC>& t)
	{
		return (1.0f - t) * v1 + t * v2;
	}

	inline vec3 slerp(const vec3& v1, const vec3& v2, float t)
	{
		quat q1;
		quat q2;
		q1 = quat::Identity();
		q2.setFromTwoVectors(v1, v2);
		return (q1.slerp(t, q2)) * v1;
	}

	inline float smoothstep(float edge0, float edge1, float x)
	{
		float t;
		t = std::clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
		return t * t * (3.0f - 2.0f * t);
	}

	enum class WrapMode
	{
		Repeat,
		Clamp,
		Natural,
	};

	enum class TickMode
	{
		Middle,
		Boundary,
	};

	template <WrapMode wrap_x, WrapMode wrap_y>
	inline void bilinear_helper(const vec2& uv, const vec2i& res, int& u0, int& u1, int& v0, int& v1, float& t0, float& t1)
	{
		static_assert(wrap_x == WrapMode::Repeat || wrap_x == WrapMode::Clamp, "Invalid wrap option.");
		static_assert(wrap_y == WrapMode::Repeat || wrap_y == WrapMode::Clamp, "Invalid wrap option.");

		vec2 uv_scaled = uv.cwiseProduct(res.cast<float>()) - vec2::Constant(0.5f);
		u0 = std::floor(uv_scaled.x());
		u1 = u0 + 1;
		v0 = std::floor(uv_scaled.y());
		v1 = v0 + 1;
		t0 = uv_scaled.x() - u0;
		t1 = uv_scaled.y() - v0;

		if constexpr (wrap_x == WrapMode::Repeat) {
			u0 = (u0 + res.x()) % res.x();
			u1 = (u1 + res.x()) % res.x();
		}
		else if (wrap_x == WrapMode::Clamp) { // Clamp
			u0 = std::clamp(u0, 0, res.x() - 1);
			u1 = std::clamp(u1, 0, res.x() - 1);
		} // else Natural

		if constexpr (wrap_y == WrapMode::Repeat) {
			v0 = (v0 + res.y()) % res.y();
			v1 = (v1 + res.y()) % res.y();
		}
		else if (wrap_y == WrapMode::Clamp) {
			v0 = std::clamp(v0, 0, res.y() - 1);
			v1 = std::clamp(v1, 0, res.y() - 1);
		} // else Natural
	}

	inline void lerp_helper(float u, int res, WrapMode wrap, TickMode tick, int& u0, int& u1, float& t)
	{
		// ASSERT(u >= 0.0f && u <= 1.0f);
		float u_scaled;
		if (tick == TickMode::Middle) {
			u_scaled = u * (float)res - 0.5f;
		}
		else {
			u_scaled = u * (float)(res - 1);
		}
		u0 = std::floor(u_scaled);
		u1 = u0 + 1;
		t = u_scaled - u0;

		if (wrap == WrapMode::Repeat) {
			u0 = (u0 + res) % res;
			u1 = (u1 + res) % res;
		}
		else if (wrap == WrapMode::Clamp) {
			u0 = std::clamp(u0, 0, res - 1);
			u1 = std::clamp(u1, 0, res - 1);
		} // else Natural
	}

	template <int N>
	inline void lerp_helper(const float* u, const int* res, WrapMode wrap, TickMode tick, int* u0, int* u1, float* t)
	{
		for (int i = 0; i < N; ++i) {
			lerp_helper(u[i], res[i], wrap, tick, u0[i], u1[i], t[i]);
		}
	}

	template <int N>
	inline void lerp_helper(const float* u, const int* res, const WrapMode* wrap, const TickMode* tick, int* u0, int* u1,
		float* t)
	{
		for (int i = 0; i < N; ++i) {
			lerp_helper(u[i], res[i], wrap[i], tick[i], u0[i], u1[i], t[i]);
		}
	}

	// Duff, Tom, et al. "Building an orthonormal basis, revisited." Journal of Computer Graphics Techniques Vol 6.1 (2017).
	inline void orthonormal_basis(const vec3& N, vec3& X, vec3& Y)
	{
		float sign = copysignf(1.0f, N.z());
		const float a = -1.0f / (sign + N.z());
		const float b = N.x() * N.y() * a;
		X = vec3(1.0f + sign * N.x() * N.x() * a, sign * b, -sign * N.x());
		Y = vec3(b, sign + N.y() * N.y() * a, -N.y());
	}

	struct Frame
	{
		Frame()
		{
			t = vec3::UnitX();
			b = vec3::UnitY();
			n = vec3::UnitZ();
		}

		Frame(const vec3& t, const vec3& b, const vec3& n) : t(t), b(b), n(n) {}

		Frame(const vec3& t, const vec3& b) : t(t), b(b) { n = t.cross(b).normalized(); }

		explicit Frame(const vec3& n) : n(n) { orthonormal_basis(n, t, b); }

		bool valid() const
		{
			return (std::abs(t.squaredNorm() - 1.0f) <= 1e-4f) && (std::abs(b.squaredNorm() - 1.0f) <= 1e-4f) &&
				(std::abs(n.squaredNorm() - 1.0f) <= 1e-4f) && (t.isOrthogonal(b, 1e-4f)) &&
				(b.isOrthogonal(n, 1e-4f)) && (n.isOrthogonal(t, 1e-4f));
		}

		vec3 to_local(const vec3& w) const { return vec3(t.dot(w), b.dot(w), n.dot(w)); }

		vec3 to_world(const vec3& w) const { return t * w.x() + b * w.y() + n * w.z(); }

		vec3 t, b, n;
	};

	inline mat4 scale_rotate_translate(const vec3& scale, const quat& rototation, const vec3& translatation)
	{
		mat4 S = mat4::Identity();
		for (int i = 0; i < 3; ++i) {
			S(i, i) = scale(i);
		}
		mat4 R = mat4::Identity();
		R.block<3, 3>(0, 0) = rototation.matrix();

		mat4 T = mat4::Identity();
		for (int i = 0; i < 3; ++i) {
			T(i, 3) = translatation(i);
		}
		return T * R * S;
	}

	inline mat4 affine_inverse(const mat4& mat)
	{
		mat4 inv;
		mat3 Ainv = mat.block(0, 0, 3, 3).inverse();
		inv.block(0, 0, 3, 3) = Ainv;
		inv.col(3).head(3) = -Ainv * mat.col(3).head(3);
		inv.row(3) = vec4(0.0f, 0.0f, 0.0f, 1.0f);
		return inv;
	}

	static inline mat4 SetPerspective(float fov, float aspect, float znear, float zfar) {
		float tan_fov_2 = tanf(fov * 0.5f);
		float yScale = 1.0 / tan_fov_2;
		float xScale = yScale / aspect;
		float zdif = znear - zfar;
		ks::mat4 Perspective;
		Perspective <<
			xScale, 0, 0, 0,
			0, yScale, 0, 0,
			0, 0, (znear + zfar) / (zdif), 2 * zfar * znear / (zdif),
			0, 0, -1, 0;
		return Perspective;
	}

	static inline mat4 SetOrtho_Camera(float l, float r, float t, float b, float n, float f) {
		ks::mat4 Ortho;
		Ortho <<
			2 / (r - l), 0, 0, -(r + l) / (r - l),
			0, 2 / (t - b), 0, -(t + b) / (t - b),
			0, 0, -1 / (f - n), n / (f - n),
			0, 0, 0, 1;
		return Ortho;
	}

	static inline mat4 SetOrtho(float l, float r, float t, float b, float n, float f) {
		ks::mat4 Ortho;
		Ortho <<
			2 / (r - l), 0, 0, -(r + l) / (r - l),
			0, 2 / (t - b), 0, -(t + b) / (t - b),
			0, 0, -2 / (f - n), -(f + n) / (f - n),
			0, 0, 0, 1;
		return Ortho;
	}

	static inline mat4 SetView(const vec3& pos, const vec3& target, const vec3& up) {
		vec3 f = target - pos;
		f = f.normalized();
		vec3 s = f.cross(up);
		s = s.normalized();
		vec3 u = s.cross(f);
		u = u.normalized();

		mat4 mat;
		mat4 m;
		// clang-format off
		m <<
			s.x(), s.y(), s.z(), 0,
			u.x(), u.y(), u.z(), 0,
			-f.x(), -f.y(), -f.z(), 0,
			0.0f, 0.0f, 0.0f, 1.0f;

		mat4 t;
		t <<
			1, 0, 0, -pos.x(),
			0, 1, 0, -pos.y(),
			0, 0, 1, -pos.z(),
			0.0f, 0.0f, 0.0f, 1.0f;
		// clang-format on
		mat = m * t;
		return mat;
	}

	inline vec3 transform_dir(const mat4& m, const vec3& d) { return m.block<3, 3>(0, 0) * d; }

	inline vec3 transform_point(const mat4& m, const vec3& p) { return (m * p.homogeneous()).head(3); }

	inline vec3 transform_normal(const mat4& m, const vec3& n)
	{
		// Consider use Transform which caches the inverse matrix.
		return (m.inverse().transpose().block<3, 3>(0, 0) * n).normalized();
	}

	inline Frame transform_frame(const mat4& m, const Frame& f)
	{
		vec3 t = transform_dir(m, f.t).normalized();
		// Need to re-project if m is not rigid...when does this happen tho
		vec3 b = transform_dir(m, f.b);
		b = (b - b.dot(t) * t).normalized();
		return Frame(t, b);
	}

	struct Transform
	{
		Transform() = default;
		explicit Transform(const mat4& m) : m(m), inv(m.inverse()) {}
		Transform(const mat4& m, const mat4& inv) : m(m), inv(inv) {}

		Transform inverse() const { return { inv, m }; }

		vec3 point(const vec3& p) const { return (m * p.homogeneous()).head(3); }

		vec4 hpoint(const vec3& p) const { return m * p.homogeneous(); }

		vec3 point_hdiv(const vec3& p) const { return (m * p.homogeneous()).hnormalized(); }

		vec3 direction(const vec3& v) const { return m.block<3, 3>(0, 0) * v; }

		vec3 normal(const vec3& n) const { return (inv.transpose().block<3, 3>(0, 0) * n).normalized(); }

		Frame frame(const Frame& f) const
		{
			vec3 t = direction(f.t).normalized();
			// Need to re-project if m is not rigid...when does this happen tho
			vec3 b = transform_dir(m, f.b);
			b = (b - b.dot(t) * t).normalized();
			return Frame(t, b);
		}

		mat4 m = mat4::Identity();
		mat4 inv = mat4::Identity();
	};

	inline Transform operator*(const Transform& a, const Transform& b)
	{
		Transform ret;
		ret.m = a.m * b.m;
		ret.inv = b.inv * a.inv;
		return ret;
	}

	inline Transform ortho_proj(float right, float left, float top, float bottom, float near_, float far_)
	{
		mat4 m = mat4::Identity();
		m(0, 0) = 2 / (right - left);
		m(1, 1) = 2 / (top - bottom);
		m(2, 2) = 2 / (near_ - far_);
		m(0, 3) = -(right + left) / (right - left);
		m(1, 3) = -(top + bottom) / (top - bottom);
		m(2, 3) = -(near_ + far_) / (near_ - far_);
		m(3, 3) = 1;
		return Transform(m);
	}

	// look at -Z axis, rotate up to Y axis, PS:use this function after normalization
	inline Transform view_proj(const vec3& pos, const vec3& lookat, const vec3& up)
	{
		vec3 right = lookat.cross(up);
		mat4 trans = mat4::Identity();
		trans(0, 3) = -pos.x();
		trans(1, 3) = -pos.y();
		trans(2, 3) = -pos.z();
		mat4 rot;
		rot << right.x(), right.y(), right.z(), 0, up.x(), up.y(), up.z(), 0, -lookat.x(), -lookat.y(), -lookat.z(), 0, 0,
			-0, 0, 1;
		return Transform(rot * trans);
	}

	inline vec3 reflect(const vec3& w, const vec3& n) { return 2.0f * n.dot(w) * n - w; }

	// eta = eta_i / eta_t
	// Make sure wi and n are on the same side.
	inline bool refract(const vec3& wi, const vec3& n, float eta, vec3& wt)
	{
		float NdotI = n.dot(wi);
		float k = 1.0f - eta * eta * (1.0f - sqr(NdotI));
		if (k < 0.0f)
			return false;
		wt = -eta * wi + (eta * NdotI - sqrt(k)) * n;
		return true;
	}

	inline mat4 SetIdentity() {
		mat4 m;
		m << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;
		return m;
	}

	/// Set a rotation matrix about the given axis by cos and sin angle theta
	inline mat4 SetRotation(const vec3& axis, float c, float s)
	{
		mat4 m;
		if (c == 1) {
			m = SetIdentity();
			return m;
		}

		float t = 1 - c;
		float tx = t * axis.x();
		float ty = t * axis.y();
		float tz = t * axis.z();
		float txy = tx * axis.y();
		float txz = tx * axis.z();
		float tyz = ty * axis.z();
		float sx = s * axis.x();
		float sy = s * axis.y();
		float sz = s * axis.z();

		m << tx * axis.x() + c, txy - sz, txz + sy, 0,
			txy + sz, ty* axis.y() + c, tyz - sx, 0,
			txz - sy, tyz + sx, tz* axis.z() + c, 0,
			0, 0, 0, 1;

		return m;
	}

	/// Set a rotation matrix about the given axis by angle theta
	inline mat4 SetRotation(const vec3& axis, float theta)
	{
		mat4 m;
		float c = (float)cos(theta);
		if (c == 1) {
			m = SetIdentity();
			return m;
		}
		float s = (float)sin(theta);
		return SetRotation(axis, c, s);
	}

	//-------------------------------------------------------------------------------

	/// Set a rotation matrix that sets <from> unit vector to <to> unit vector
	inline mat4 SetRotation(const vec3& from, const vec3& to)
	{
		float c = from.dot(to);
		if (c > 0.999999) {
			mat4 m = SetIdentity();
			return m;
		}
		float s = (float)sqrt(1 - c * c);
		vec3 axis = from.cross(to);
		return SetRotation(axis, c, s);
	}

	/// Sets the translation component of the matrix
	inline mat4 SetTrans(const vec3& move)
	{
		mat4 m;
		m << 1, 0, 0, move.x(),
			0, 1, 0, move.y(),
			0, 0, 1, move.z(),
			0, 0, 0, 1;
		return m;
	}

	static mat4 SetRotationX(float theta) { return SetRotation(vec3(1, 0, 0), theta); }
	/// Set as rotation matrix around y axis in radians
	static mat4 SetRotationY(float theta) { return SetRotation(vec3(0, 1, 0), theta); }
	/// Set as rotation matrix around z axis in radians
	static mat4 SetRotationZ(float theta) { return SetRotation(vec3(0, 0, 1), theta); }
	/// Returns a rotation matrix around x axis in radians
	static mat4 MatrixRotationX(float theta) { return SetRotationX(theta); }
	/// Returns a rotation matrix around y axis in radians
	static mat4 MatrixRotationY(float theta) { return SetRotationY(theta); }
	/// Returns a rotation matrix around z axis in radians
	static mat4 MatrixRotationZ(float theta) { return SetRotationZ(theta); }
	/// Returns a translation matrix with no rotation or scale
	static mat4 MatrixTrans(const vec3& move) {
		mat4 m = SetTrans(move);
		return m;
	}

} // namespace ks