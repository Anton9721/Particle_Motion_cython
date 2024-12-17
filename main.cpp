#include "Particle.hpp"
#include "Vector3D.hpp"
#include "Writer.hpp"
#include "Pusher.hpp"
#include "RadiationSpectrum.hpp"

#include <omp.h>

int main()
{
    double pi = std::numbers::pi;

    std::string file_trajectory = "data_particles_crossEMField.txt";
    std::string file_spectrum = "radiation_spectrum_crossEMField.txt";
    std::string file_color_map = "color_map_crossEMField.txt";

    //Задание параметров частицы и полей
    std::vector<Particle>
        particles = {
            {Vector3D(0.0, 0.0, 0.0), Vector3D(0.1, 0.1, 0.0), 1.0, -1.0}, // position, velocity (<1), mass, charge
        };

    Vector3D electric_field(0.0, 0.0, 0.0);
    Vector3D magnetic_field(0.0, 1.0, 0.0);

    // double omega_f = 1.0;
    // double e_initial_phase = 0.0;
    // double b_initial_phase = 0.0;
    // OscillatingEMField em_field(electric_field, magnetic_field, omega_f, e_initial_phase, b_initial_phase);


    // Нахождение траектории и скоростей
    PusherBorisRR solver;
    WriterTXT writer_trajectory(file_trajectory);

    double time = 0.0;
    double time_end = 170 * 2 * pi;
    double dt = 0.01;
    int steps = int((time_end - time) / dt);

    GaussEMField em_field(electric_field, magnetic_field, time, time_end);

    for (int i = 0; i < steps; ++i)
    {
        solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
        writer_trajectory.write(std::span<Particle>(particles.data(), particles.size()), time);
        time += dt;
    }

    //Поиск спектра
    WriterSpectrum writer_spectrum(file_spectrum);

    double phi = 0;
    double theta = pi/2;
    double w_start = 0.1;
    double w_end = 4.0;
    double dw = 0.05;
    int steps_w = static_cast<int>((w_end - w_start) / dw) + 1;

    double dtheta = 0.05;
    int steps_theta = static_cast<int>(2 * pi / dtheta) + 1;

    TrajectoryReader reader(file_trajectory);
    std::span<const Vector3D> trajectory_s = reader.getTrajectory();
    std::span<const Vector3D> velocity_s = reader.getVelocity();
    std::span<const double> time_s = reader.getTime();
    std::vector<double> spectrum_data;

    RadiationSpectrum spectrum(trajectory_s, velocity_s, time_s);
    spectrum_data = spectrum.calculateSpectrum(phi, theta, w_start, w_end, dw);
    // writer_spectrum.write(spectrum_data, w_start, w_end, dw);


    // Цветовая карта
    WriterSpectrumMap writer_map(file_color_map);
    std::vector<std::vector<double>> color_map(steps_theta, std::vector<double>(steps_w));
    for (int i = 0; i < steps_theta; ++i)
    {
        double theta_map = i * 2 * pi / (steps_theta - 1);
        for (int j = 0; j < steps_w; ++j)
        {
            double w_map = j * w_end / (steps_w - 1);
            color_map[i][j] = spectrum.calculateIntensity(phi, theta_map, w_map);
        }
    }

    writer_map.write(color_map);

    return 0;
}

// mkdir build
// cd build
// cmake ..
// cmake --build .
// cd..
// python setup.py build_ext-- inplace

//  .\Debug\ParticleMotion.exe