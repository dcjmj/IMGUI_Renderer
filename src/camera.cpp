#include "camera.h"
namespace ks
{

    static inline mat4 look_at(const vec3& position, const vec3& target, vec3 up)
    {
        vec3 back = (position - target).normalized();
        vec3 right = up.cross(back);

        if (right.squaredNorm() < 1e-8f) {
            orthonormal_basis(back, right, up);
        }
        else {
            right.normalize();
            up = back.cross(right);
        }
        mat4 mat;
        // clang-format off
        mat <<
            right.x(), up.x(), back.x(), position.x(),
            right.y(), up.y(), back.y(), position.y(),
            right.z(), up.z(), back.z(), position.z(),
            0.0f, 0.0f, 0.0f, 1.0f;
        // clang-format on
        return mat;
    }

    static inline mat4 look_at_view(const vec3& position, const vec3& target, vec3 up)
    {
        vec3 back = (position - target).normalized();
        vec3 right = up.cross(back);

        if (right.squaredNorm() < 1e-8f) {
            orthonormal_basis(back, right, up);
        }
        else {
            right.normalize();
            up = back.cross(right);
        }
        mat4 mat;
        // clang-format off
        mat <<
            right.x(), right.y(), right.z(), -right.dot(position),
            up.x(), up.y(), up.z(), -up.dot(position),
            back.x(), back.y(), back.z(), -back.dot(position),
            0.0f, 0.0f, 0.0f, 1.0f;
        // clang-format on
        return mat;
    }

    static inline mat4 rev_inf_projection(float vfov, float aspect, float near_clip = 0.001f)
    {
        // Reverse inf projection.
        float cot_half_vfov = 1.0f / std::tan(vfov * 0.5f);
        // Vulkan NDC space y axis pointing downward (same as screen space).
        mat4 proj;
        float far = 100.f;
        float near = 0.1f;
        // clang-format off
        proj <<
            cot_half_vfov / aspect, 0.0f, 0.0f, 0.0f,
            0.0f, -cot_half_vfov, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, near,
            0.0f, 0.0f, -1.0f, 0.0f;
        // clang-format on
        return proj;
    }

    static inline mat4 rev_orthographic(float left, float right, float bottom, float top, float near, float far)
    {
        mat4 proj;
        // clang-format off
        proj <<
            2.0f / (right - left), 0.0f, 0.0f, -(right + left) / (right - left),
            0.0f, 2.0f / (bottom - top), 0.0f, -(top + bottom) / (bottom - top),
            0.0f, 0.0f, 1.0f / (far - near), far / (far - near),
            0.0f, 0.0f, 0.0f, 1.0f;
        // clang-format on
        return proj;
    }

} // namespace ks