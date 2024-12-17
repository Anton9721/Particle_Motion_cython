#pragma once

#include <omp.h>
#include <vector>
#include <span>
#include <cmath>
#include <iostream>
#include <complex>
#include <fstream>
#include "Vector3D.hpp" 
#include "Particle.hpp"

class RadiationSpectrum
{
public:
    RadiationSpectrum(std::span<const Vector3D> &trajectory, std::span<const Vector3D> &velocity, std::span<const double> &time)
        : trajectory_(trajectory), velocity_(velocity), time_(time) {}

    double calculateIntensity(double phi, double theta, double omega) const;
    std::vector<double> calculateSpectrum(double phi, double theta, double w_start, double w_end, double dw);
    std::vector<double> calculateSpectrumMap(double phi, double theta, double w_start, double w_end, double dw, double dtheta);

private:
    std::span<const Vector3D> trajectory_;
    std::span<const Vector3D> velocity_;
    std::span<const double> time_;

    double ksi(int j, const Vector3D &normal) const
    {
        return time_[j] - dot(normal, trajectory_[j]);
    }

    Vector3D normal(double phi, double theta) const
    {
        return Vector3D(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    }
};
