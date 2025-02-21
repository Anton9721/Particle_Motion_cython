#include "Vector3D.hpp"
#include <cmath>

Vector3D operator+(const Vector3D &v1, const Vector3D &v2)
{
    return Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector3D operator-(const Vector3D &v1, const Vector3D &v2)
{
    return Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

Vector3D operator*(const Vector3D &v, double scalar)
{
    return Vector3D(v.x * scalar, v.y * scalar, v.z * scalar);
}

Vector3D operator*(double scalar, const Vector3D &v)
{
    return Vector3D(v.x * scalar, v.y * scalar, v.z * scalar);
}

Vector3D operator/(const Vector3D &v, double scalar)
{
    return scalar != 0 ? Vector3D(v.x / scalar, v.y / scalar, v.z / scalar) : v;
}

Vector3D &operator+=(Vector3D &v1, const Vector3D &v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

Vector3D &operator-=(Vector3D &v1, const Vector3D &v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

Vector3D &operator*=(Vector3D &v, double scalar)
{
    v.x *= scalar;
    v.y *= scalar;
    v.z *= scalar;
    return v;
}

Vector3D &operator/=(Vector3D &v, double scalar)
{
    if (scalar != 0)
    {
        v.x /= scalar;
        v.y /= scalar;
        v.z /= scalar;
    }
    return v;
}

Vector3D cross(const Vector3D &v1, const Vector3D &v2)
{
    return Vector3D(
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

double dot(const Vector3D &v1, const Vector3D &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double abs(const Vector3D &v)
{
    return std::sqrt(dot(v, v));
}

std::ostream &operator<<(std::ostream &os, const Vector3D &v)
{
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}
