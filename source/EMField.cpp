#include "EMField.hpp"

Vector3D CrossEMField::get_electric_field(const Vector3D &r, double &t) const
{
    return electric_field_;
}

Vector3D CrossEMField::get_magnetic_field(const Vector3D &r, double &t) const
{
    return magnetic_field_;
}

Vector3D OscillatingEMField::get_electric_field(const Vector3D &r, double &t) const
{
    return electric_field_ * std::sin(omega * t + initial_phase_e);
}

Vector3D OscillatingEMField::get_magnetic_field(const Vector3D &r, double &t) const
{
    return magnetic_field_ * std::sin(omega * t + initial_phase_b);
}

Vector3D GaussEMField::get_electric_field(const Vector3D &r, double &t) const
{
    double gauss_factor = std::exp(-std::pow(t - mu_, 2) / (2 * sigma_ * sigma_));
    return electric_field_;
}

Vector3D GaussEMField::get_magnetic_field(const Vector3D &r, double &t) const
{
    double gauss_factor = std::exp(-std::pow(t - mu_, 2) / (2 * sigma_ * sigma_));
    return magnetic_field_;
}
