#pragma once

struct Point3D
{
    double x;
    double y;
    double z;

    Point3D() = default;
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

