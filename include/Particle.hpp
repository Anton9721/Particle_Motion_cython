#pragma once
#include "Vector3D.hpp"

struct Particle
{
    Vector3D position;
    Vector3D velocity;
    double mass;
    double charge;

    Particle(Vector3D pos, Vector3D vel, double m, double q)
        : position(pos), velocity(vel), mass(m), charge(q) {}
};
