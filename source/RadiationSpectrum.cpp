#include "RadiationSpectrum.hpp"
#include <cmath>
#include <fstream>
#include <complex>
#include <omp.h>
#include <vector>

double RadiationSpectrum::calculateIntensity(double phi, double theta, double omega) const
{
    Vector3D real_intensity(0.0, 0.0, 0.0);
    Vector3D imag_intensity(0.0, 0.0, 0.0);
    Vector3D normal_ = normal(phi, theta);

    for (int i = 0; i < time_.size() - 1; i++)
    {
        double delta_t = time_[i + 1] - time_[i];
        double delta_ksi = ksi(i + 1, normal_) - ksi(i, normal_);
        double average_ksi = (ksi(i + 1, normal_) + ksi(i, normal_)) / 2.0;
        Vector3D delta_v = velocity_[i + 1] - velocity_[i];
        Vector3D average_v = (velocity_[i + 1] + velocity_[i]) / 2.0;

        double phase_real = omega * average_ksi;
        double phase_imag = omega * delta_ksi / 2.0;

        Vector3D delta_real = 2 * average_v * std::sin(phase_imag);
        Vector3D delta_imag = delta_v * (std::sin(phase_imag) / phase_imag - std::cos(phase_imag));

        delta_real = delta_real * std::cos(phase_real) - delta_imag * std::sin(phase_real);
        delta_imag = delta_real * std::sin(phase_real) + delta_imag * std::cos(phase_real);

        real_intensity += delta_real * (delta_t / delta_ksi);
        imag_intensity += delta_imag * (delta_t / delta_ksi);
    }

    double J1_r = -(real_intensity.x * std::cos(phi) + real_intensity.y * std::sin(phi)) * std::cos(theta) -
                  real_intensity.z * std::sin(theta);
    double J1_im = -(imag_intensity.x * std::cos(phi) + imag_intensity.y * std::sin(phi)) * std::cos(theta) -
                   imag_intensity.z * std::sin(theta);
    double J2_r = real_intensity.x * std::sin(phi) - real_intensity.y * std::cos(phi);
    double J2_im = imag_intensity.x * std::sin(phi) - imag_intensity.y * std::cos(phi);

    double J1 = J1_r * J1_r + J1_im * J1_im;
    double J2 = J2_r * J2_r + J2_im * J2_im;
    double J = J1 + J2;

    return J;
}

std::vector<double> RadiationSpectrum::calculateSpectrum(double phi, double theta, double w_start, double w_end, double dw)
{
    int num_steps = static_cast<int>((w_end - w_start) / dw) + 1;
    std::vector<double> results(num_steps);

    #pragma omp parallel for
    for (int i = 0; i < num_steps; ++i)
    {
        double w = w_start + i * dw;
        results[i] = calculateIntensity(phi, theta, w);
    }   

    return results;
}

