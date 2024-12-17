#pragma once
#include "Particle.hpp"
#include "Vector3D.hpp"
#include "cmath"

class EMField
{
public:
    virtual Vector3D get_electric_field(const Vector3D &r, double &t) const = 0;
    virtual Vector3D get_magnetic_field(const Vector3D &r, double &t) const = 0;

    virtual ~EMField() = default;
};

// Создание скрещенных электромагнитных полей
class CrossEMField : public EMField
{
public:
    Vector3D get_electric_field(const Vector3D &r, double &t) const override;
    Vector3D get_magnetic_field(const Vector3D &r, double &t) const override;

    CrossEMField(Vector3D &electric_field, Vector3D &magnetic_field)
        : electric_field_(electric_field), magnetic_field_(magnetic_field)
    {
    }

private:
    Vector3D electric_field_;
    Vector3D magnetic_field_;
};

// Создание переменных (осциллирующих) электромагнитных полей
class OscillatingEMField : public EMField
{
public:
    Vector3D get_electric_field(const Vector3D &r, double &t) const override;
    Vector3D get_magnetic_field(const Vector3D &r, double &t) const override;

    OscillatingEMField(Vector3D &electric_field, Vector3D &magnetic_field, double &omega, double &initial_phase_e, double &initial_phase_b)
        : electric_field_(electric_field), magnetic_field_(magnetic_field), omega(omega), initial_phase_e(initial_phase_e), initial_phase_b(initial_phase_b)
    {
    }

private:
    Vector3D electric_field_;
    Vector3D magnetic_field_;
    double omega;
    double initial_phase_e;
    double initial_phase_b;
};

// Создание скрещенных выключающихся электромагнитных полей
class GaussEMField : public EMField
{
public:
    Vector3D get_electric_field(const Vector3D &r, double &t) const override;
    Vector3D get_magnetic_field(const Vector3D &r, double &t) const override;

    GaussEMField(Vector3D &electric_field, Vector3D &magnetic_field, double &time_start, double &time_end)
        : electric_field_(electric_field), magnetic_field_(magnetic_field), time_start_(time_start), time_end_(time_end)
    {
    }

private:
    Vector3D electric_field_;
    Vector3D magnetic_field_;
    double time_start_;
    double time_end_;
    const double mu_ = (time_start_ + time_end_) / 2.0;
    const double sigma_ = (time_end_ - time_start_) / 6.0;
};