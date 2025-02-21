#pragma once

#include "Particle.hpp"
#include "Vector3D.hpp"
#include "Writer.hpp"
#include "Pusher.hpp"
#include "RadiationSpectrum.hpp"

#include <vector>
#include <string>
#include <span>
#include <omp.h>

std::vector<double> find_trajectory(std::string solver_name, std::string field_name,
                                    std::vector<double> position, std::vector<double> velocity, double mass, double charge,
                                    std::vector<double> electric_field, std::vector<double> magnetic_field,
                                    double time_start, double time_end, double dt, double omega_f, double e_initial_phase, double b_initial_phase);

std::vector<double> spectrum(std::vector<double> time, std::vector<double> trajectory, std::vector<double> velocity,
                             double phi, double theta, double w_start, double w_end, double dw);

std::vector<std::vector<double>> color_map(std::vector<double> time, std::vector<double> trajectory, std::vector<double> velocity,
                                           double phi, double theta, double w_start, double w_end, double dw, double dtheta);
