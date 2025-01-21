#pragma once
#include "maths.h"
#include <memory>
#include <iostream>
#include <cstdio>

namespace ks
{
	class CameraBase
	{
	private:
		mat4 matrix, matrixInverse;
		mat3 normal;
	protected:
		mat4 proj;
		mat4 view;
	public:
		const mat4& GetMatrix() const { return matrix; }
		const mat4& GetViewMatrix() const { return view; }
		const mat4& GetProjMatrix() const { return proj; }
		const mat3& GetNormalMatrix() const { return normal; }
		const mat4& GetMatrixInverse() const { return matrixInverse; }

		void SetMatrix(const mat4& viewMatrix, const mat4& projMatrix) { view = viewMatrix; proj = projMatrix; }
		void UpdateViewMatrix() { UpdateView(); UpdateMatrix(); }
	private:
		void UpdateProjection() {
			ComputeProjectionMatrix();
		}
		void UpdateView() {
			ComputeViewMatrix();
			mat3 view3;
			view3 = view.block<3, 3>(0, 0);
			normal = view3.inverse().transpose();
			UpdateMatrix();
		}
		void UpdateMatrix() {
			matrix = proj * view;
			matrixInverse = matrix.inverse();
		}
	protected:
		virtual void ComputeViewMatrix() = 0;
		virtual void ComputeProjectionMatrix() = 0;
		void UpdateProjectionMatrix() { UpdateProjection(); UpdateMatrix(); }

		void UpdateAllMatrices() { UpdateView(); UpdateProjection(); UpdateMatrix(); }
	};

	//-------------------------------------------------------------------------------

	class ProjectionCamera : public CameraBase
	{
	protected:
		float fov;
		float aspect;
		float znear, zfar;
		
	public:
		float l, r, b, t, z_near, z_far;
		ProjectionCamera() : fov(20), aspect(1.5), znear(0.02f), zfar(150.0f), l(-12), r(12), b(-8), t(8), z_near(0), z_far(512.f) {}

		void SetAspect(float asp) { aspect = asp; UpdateProjectionMatrix(); }

	protected:
		void ComputeProjectionMatrix() {
			//proj = SetPerspective(DEG2RAD(fov), aspect, znear, zfar);
			proj = SetOrtho_Camera(l, r, t, b, z_near,z_far);
		}
	};

	//-------------------------------------------------------------------------------

	class Camera : public ProjectionCamera
	{
	protected:
		float distance;

		vec3 target;
		float height;
	public:
		vec2 rot;
		//Camera() : distance(15.0f), rot(0, 0), target(0, 0, 0), height(0) { UpdateAllMatrices(); }
		//Camera() : distance(5.0f), rot(0, 0), target(0, 0, 0), height(0) { UpdateAllMatrices(); }
		//Camera() : distance(1.5f), rot(0,0), target(0,0,0), height(0) { UpdateAllMatrices(); }	// regular
		//Camera() : distance(5.f), rot(0, 0), target(0, 0, 0), height(0) { UpdateAllMatrices(); }
		//Camera() : distance(16.0f), rot(3.14f * 0.2f, 3.14f * 0.12f), target(1.2, -0.7, -2.6), height(0) { UpdateAllMatrices(); } // sweater
		Camera() : distance(20.0f), rot(0, 0), target(0, 0, 0), height(0) { UpdateAllMatrices(); }

		void SetTarget(const vec3& p) { target = p;	UpdateViewMatrix(); }
		void AddRotation(float x, float y) { rot.x() += x; rot.y() += y;	UpdateViewMatrix(); }
		void AddDistance(float d) { distance += d; distance = fmax(distance, 0.1f);	l += d*1.5; r -= d*1.5; t -= d; b += d;		UpdateViewMatrix(); }
		void SetDistance(float d) { distance = d; distance = fmax(distance, 0.1f);			UpdateViewMatrix(); }
		float GetDistance() const { return distance; }
		void AddHeight(float h) { height += h;			UpdateViewMatrix(); }
		ks::vec3 GetCameraPosition() {
			// Invert the view matrix to get the camera's world transform
			ks::mat4 inverseView = view.inverse();

			// The camera position is in the last column of the inverse view matrix
			ks::vec3 cameraPosition = inverseView.block<3, 1>(0, 3); // Extracts the translation part
			//std::cout << cameraPosition << std::endl;
			return cameraPosition;
		}


	protected:
		virtual void ComputeViewMatrix() {
			view = MatrixTrans(vec3(0, 0, -distance)) * MatrixRotationX(rot.y()) * MatrixRotationY(rot.x()) * MatrixTrans(-target - vec3(0, height, 0));
		}
	};

} // namespace ks