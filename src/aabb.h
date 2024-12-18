#pragma once

#include "maths.h"
namespace ks
{
    template <typename T>
    inline T lerp_mine(T A, T B, float t){
        return float(A) + t * (float(B) - float(A));
    } 
    
    struct AABB2
    {
        AABB2() = default;
        AABB2(const vec2& pMin_, const vec2& pMax_) : pMin(pMin), pMax(pMax) {}

        void expand(const vec2& point)
        {
            pMin.x() = (std::min)(pMin.x(), point.x());
            pMin.y() = (std::min)(pMin.y(), point.y());
            pMax.x() = (std::max)(pMax.x(), point.x());
            pMax.y() = (std::max)(pMax.y(), point.y());
        }

        void expand(const AABB2& aabb)
        {
            pMin.x() = (std::min)(pMin.x(), aabb.pMin.x());
            pMin.y() = (std::min)(pMin.y(), aabb.pMin.y());

            pMax.x() = (std::max)(pMax.x(), aabb.pMax.x());
            pMax.y() = (std::max)(pMax.y(), aabb.pMax.y());
        }

        bool contain(const vec2& point) const
        {
            return point.x() >= pMin.x() && point.x() <= pMax.x() && point.y() >= pMin.y() && point.y() <= pMax.y();
        }

        bool contain(const AABB2& bound) const
        {
            return bound.pMin.x() >= pMin.x() && bound.pMax.x() <= pMax.x() && bound.pMin.y() >= pMin.y() &&
                bound.pMax.y() <= pMax.y();
        }

        bool isEmpty() const { return pMin.x() > pMax.x() || pMin.y() > pMax.y(); }

        vec2 center() const { return 0.5f * (pMin + pMax); }

        vec2 extents() const
        {
            if (isEmpty()) {
                return vec2::Zero();
            }
            return pMax - pMin;
        }

        float offset(const vec2& point, uint32_t dim) const
        {
            float ext = extents()[dim];
            if (ext == 0.0) {
                return 0.0;
            }
            return (point[dim] - pMin[dim]) / ext;
        }

        vec2 offset(const vec2& point) const
        {
            vec2 ext = extents();
            vec2 o = (point - pMin).cwiseQuotient(ext);
            for (int i = 0; i < 2; ++i) {
                if (ext[i] == 0.0f)
                    o[i] = 0.0f;
            }
            return o;
        }

        vec2 lerp(const vec2& t) const
        {
            return vec2(lerp_mine(pMin.x(), pMax.x(), t.x()), lerp_mine(pMin.y(), pMax.y(), t.y()));
        }

        uint32_t largestAxis() const
        {
            vec2 exts = extents();
            if (exts.x() >= exts.y())
                return 0;
            else
                return 1;
        }

        float area() const
        {
            vec2 exts = extents();
            return exts.x() * exts.y();
        }

        const vec2& operator[](uint32_t index) const
        {
            return (reinterpret_cast<const vec2*>(this))[index];
        }

        vec2& operator[](uint32_t index)
        {
            return (reinterpret_cast<vec2*>(this))[index];
        }

        vec2 corner(int i) const
        {
            const vec2* c = (const vec2*)this;
            return vec2(c[i & 1].x(), c[(i >> 1) & 1].y());
        }

        vec2 pMin = vec2::Constant(inf);
        vec2 pMax = vec2::Constant(-inf);
    };

    inline bool intersectBool(const AABB2& b0, const AABB2& b1)
    {
        return !(b0.pMin.x() > b1.pMax.x() || b1.pMin.x() > b0.pMax.x() || b0.pMin.y() > b1.pMax.y() || b1.pMin.y() > b0.pMax.y());
    }

    inline AABB2 intersect(const AABB2& b0, const AABB2& b1)
    {
        AABB2 ret;
        ret.pMin = b0.pMin.cwiseMax(b1.pMin);
        ret.pMax = b0.pMax.cwiseMin(b1.pMax);
        if (ret.pMin.x() > ret.pMax.x() || ret.pMin.y() > ret.pMax.y()) {
            return AABB2();
        }
        return ret;
    }

    inline AABB2 join(const AABB2& b0, const AABB2& b1)
    {
        AABB2 ret;
        ret.pMin = b0.pMin.cwiseMin(b1.pMin);
        ret.pMax = b0.pMax.cwiseMax(b1.pMax);
        return ret;
    }

    struct AABB3
    {
        AABB3() = default;

        void reset() {
            float float_min = (std::numeric_limits<float>::min)(),
                float_max = (std::numeric_limits<float>::max)();
            pMin = vec3(float_max, float_max, float_max);
            pMax = vec3(float_min, float_min, float_min);
        }
        
        AABB3(const vec3& pMin, const vec3& pMax) : pMin(pMin), pMax(pMax) {}

        void expand(const vec3& point)
        {
            pMin.x() = (std::min)(pMin.x(), point.x());
            pMin.y() = (std::min)(pMin.y(), point.y());
            pMin.z() = (std::min)(pMin.z(), point.z());

            pMax.x() = (std::max)(pMax.x(), point.x());
            pMax.y() = (std::max)(pMax.y(), point.y());
            pMax.z() = (std::max)(pMax.z(), point.z());
        }

        void expand(const AABB3& aabb)
        {
            pMin.x() = (std::min)(pMin.x(), aabb.pMin.x());
            pMin.y() = (std::min)(pMin.y(), aabb.pMin.y());
            pMin.z() = (std::min)(pMin.z(), aabb.pMin.z());

            pMax.x() = (std::max)(pMax.x(), aabb.pMax.x());
            pMax.y() = (std::max)(pMax.y(), aabb.pMax.y());
            pMax.z() = (std::max)(pMax.z(), aabb.pMax.z());
        }

        bool isEmpty() const { return pMin.x() > pMax.x() || pMin.y() > pMax.y() || pMin.z() > pMax.z(); }

        vec3 center() const { return 0.5f * (pMin + pMax); }

        vec3 extents() const { return pMax - pMin; }

        float offset(const vec3& point, uint32_t dim) const
        {
            float ext = extents()[dim];
            if (ext == 0.0f) {
                return 0.0f;
            }
            return (point[dim] - pMin[dim]) / ext;
        }

        vec3 offset(const vec3& point) const
        {
            vec3 ext = extents();
            vec3 o = (point - pMin).cwiseQuotient(ext);
            for (int i = 0; i < 3; ++i) {
                if (ext[i] == 0.0f)
                    o[i] = 0.0f;
            }
            return o;
        }

        vec3 corner(int i) const
        {
            const vec3* c = (const vec3*)this;
            return vec3(c[i & 1].x(), c[(i & 2) >> 1].y(), c[(i & 4) >> 2].z());
        }

        vec3 lerp(const vec3& t) const
        {
            return vec3(lerp_mine(pMin.x(), pMax.x(), t.x()), lerp_mine(pMin.y(), pMax.y(), t.y()),
                lerp_mine(pMin.z(), pMax.z(), t.z()));
        }

        uint32_t largestAxis() const
        {
            vec3 exts = extents();
            if (exts.x() >= exts.y() && exts.x() >= exts.z())
                return 0;
            else if (exts.y() >= exts.x() && exts.y() >= exts.z())
                return 1;
            else
                return 2;
        }

        float surfaceArea() const
        {
            vec3 exts = extents();
            return 2.0f * (exts.x() * exts.y() + exts.y() * exts.z() + exts.x() * exts.z());
        }

        float volume() const
        {
            vec3 exts = extents();
            return exts.x() * exts.y() * exts.z();
        }

        bool contain(const vec3& p) const
        {
            return pMin.x() <= p.x() && pMin.y() <= p.y() && pMin.z() <= p.z() && pMax.x() >= p.x() && pMax.y() >= p.y() &&
                pMax.z() >= p.z();
        }

        bool contain(const AABB3& other) const
        {
            return pMin.x() <= other.pMin.x() && pMin.y() <= other.pMin.y() && pMin.z() <= other.pMin.z() &&
                pMax.x() >= other.pMax.x() && pMax.y() >= other.pMax.y() && pMax.z() >= other.pMax.z();
        }

        const vec3& operator[](uint32_t index) const
        {
            return (reinterpret_cast<const vec3*>(this))[index];
        }

        vec3& operator[](uint32_t index)
        {
            return (reinterpret_cast<vec3*>(this))[index];
        }

        vec3 pMin = vec3::Constant(inf);
        vec3 pMax = vec3::Constant(-inf);
    };

    inline bool intersectBool(const AABB3& b0, const AABB3& b1)
    {
        return !(b0.pMin.x() > b1.pMax.x() || b1.pMin.x() > b0.pMax.x() || b0.pMin.y() > b1.pMax.y() || b1.pMin.y() > b0.pMax.y() ||
            b0.pMin.x() > b1.pMax.z() || b1.pMin.z() > b0.pMax.z());
    }

    inline AABB3 intersect(const AABB3& b0, const AABB3& b1)
    {
        AABB3 ret;
        ret.pMin = b0.pMin.cwiseMax(b1.pMin);
        ret.pMax = b0.pMax.cwiseMin(b1.pMax);
        if (ret.pMin.x() > ret.pMax.x() || ret.pMin.y() > ret.pMax.y() || ret.pMin.z() > ret.pMax.z()) {
            return AABB3();
        }
        return ret;
    }

    inline AABB3 join(const AABB3& b0, const AABB3& b1) { return AABB3(b0.pMin.cwiseMin(b1.pMin), b0.pMax.cwiseMax(b1.pMax)); }

    inline AABB3 transform_aabb(const Transform& t, const AABB3& b)
    {
        const vec3* c = (const vec3*)&b;
        AABB3 out;
        for (int i = 0; i < 8; ++i) {
            vec3 corner(c[i & 1].x(), c[(i & 2) >> 1].y(), c[(i & 4) >> 2].z());
            out.expand(t.point(corner));
        }
        return out;
    }

} // namespace ks