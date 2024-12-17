#pragma once

#include <span>
#include <cmath>
#include <numbers>
#include "Particle.hpp"
#include "Vector3D.hpp"
#include "EMField.hpp"


class Pusher
{
public:
    virtual void push(std::span<Particle> particles, const EMField &field, double t, double dt) const = 0;
    virtual ~Pusher() = default;                                                                          
};

// Решение методом Эйлера
class PusherEuler : public Pusher
{
public:
    void push(std::span<Particle> particles, const EMField &field, double t, double dt) const override;
};

// Решение методом Рунге-Кутты 4 порядка
class RungeKutta4 : public Pusher
{
public:
    void push(std::span<Particle> particles, const EMField &field, double t, double dt) const override;
    Vector3D get_Lorentz_force(const Particle &particle, const EMField &field, const Vector3D &r, const Vector3D &velocity, double t) const;
};

// Решение методом Бориса без учета радиационной силы
class PusherBoris : public Pusher
{
public:
    void push(std::span<Particle> particles, const EMField &field, double t, double dt) const override;
};

// Решение методом Бориса с учетом радиационной силы
class PusherBorisRR : public Pusher
{
public:
    void push(std::span<Particle> particles, const EMField &field, double t, double dt) const override;
};

