#pragma once

#include "Point3D.hpp"
#include <iostream>

typedef Point3D Vector3D;

Vector3D operator+(const Vector3D &v1, const Vector3D &v2);
Vector3D operator-(const Vector3D &v1, const Vector3D &v2);
Vector3D operator*(const Vector3D &v, double scalar);
Vector3D operator*(double scalar, const Vector3D &v);
Vector3D operator/(const Vector3D &v, double scalar);

Vector3D &operator+=(Vector3D &v1, const Vector3D &v2);
Vector3D &operator-=(Vector3D &v1, const Vector3D &v2);
Vector3D &operator*=(Vector3D &v, double scalar);
Vector3D &operator/=(Vector3D &v, double scalar);

Vector3D cross(const Vector3D &v1, const Vector3D &v2);
double dot(const Vector3D &v1, const Vector3D &v2);
double abs(const Vector3D &v);

std::ostream &operator<<(std::ostream &os, const Vector3D &vec);

