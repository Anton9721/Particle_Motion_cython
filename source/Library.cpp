# include "Library.hpp"

std::vector<double> find_trajectory(std::string solver_name, std::string field_name,
                                                std::vector<double> position, std::vector<double> velocity, double mass, double charge,
                                                std::vector<double> electric_field, std::vector<double> magnetic_field,
                                                double time_start, double time_end, double dt, double omega_f, double e_initial_phase, double b_initial_phase)
{
    Vector3D pos(position[0], position[1], position[2]);
    Vector3D vel(velocity[0], velocity[1], velocity[2]);
    Vector3D e_field(electric_field[0], electric_field[1], electric_field[2]);
    Vector3D b_field(magnetic_field[0], magnetic_field[1], magnetic_field[2]);

    // double omega_f = 1.0;
    // double e_initial_phase = 0.0;
    // double b_initial_phase = 0.0;

    int steps = int((time_end - time_start) / dt);
    double time = time_start;
    std::vector<double> data_result;

    std::vector<Particle>
        particles = {
            {pos, vel, mass, charge}, // position, velocity (<1), mass, charge
        };

    WriterPyhton writer_trajectory;

    // Перекрестное поле
    if (solver_name == "PusherEuler" && field_name == "CrossEMField")
    {
        PusherEuler solver;
        CrossEMField em_field(e_field, b_field);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "RungeKutta4" && field_name == "CrossEMField")
    {
        RungeKutta4 solver;
        CrossEMField em_field(e_field, b_field);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }

    else if (solver_name == "PusherBoris" && field_name == "CrossEMField")
    {
        PusherBoris solver;
        CrossEMField em_field(e_field, b_field);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "PusherBorisRR" && field_name == "CrossEMField")
    {
        PusherBorisRR solver;
        CrossEMField em_field(e_field, b_field);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }

    // Поле по Гауссу
    else if (solver_name == "PusherEuler" && field_name == "GaussEMField")
    {
        PusherEuler solver;
        GaussEMField em_field(e_field, b_field, time_start, time_end);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "RungeKutta4" && field_name == "GaussEMField")
    {
        RungeKutta4 solver;
        GaussEMField em_field(e_field, b_field, time_start, time_end);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "PusherBoris" && field_name == "GaussEMField")
    {
        PusherBoris solver;
        GaussEMField em_field(e_field, b_field, time_start, time_end);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "PusherBorisRR" && field_name == "GaussEMField")
    {
        PusherBorisRR solver;
        GaussEMField em_field(e_field, b_field, time_start, time_end);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }

    // Осциллирующее поле
    else if (solver_name == "PusherEuler" && field_name == "OscillatingEMField")
    {
        PusherEuler solver;
        OscillatingEMField em_field(e_field, b_field, omega_f, e_initial_phase, b_initial_phase);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }

    else if (solver_name == "RungeKutta4" && field_name == "OscillatingEMField")
    {
        RungeKutta4 solver;
        OscillatingEMField em_field(e_field, b_field, omega_f, e_initial_phase, b_initial_phase);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "PusherBoris" && field_name == "OscillatingEMField")
    {
        PusherBoris solver;
        OscillatingEMField em_field(e_field, b_field, omega_f, e_initial_phase, b_initial_phase);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    else if (solver_name == "PusherBorisRR" && field_name == "OscillatingEMField")
    {
        PusherBorisRR solver;
        OscillatingEMField em_field(e_field, b_field, omega_f, e_initial_phase, b_initial_phase);
        for (int i = 0; i < steps; ++i)
        {
            solver.push(std::span<Particle>(particles.data(), particles.size()), em_field, time, dt);
            std::vector<double> res = writer_trajectory.writer_numpy(std::span<Particle>(particles.data(), particles.size()), time);
            for (const auto &member : res)
            {
                data_result.push_back(member);
            }

            time += dt;
        }
    }
    return data_result;
}

std::vector<double> spectrum(std::vector<double> time, std::vector<double> trajectory, std::vector<double> velocity, 
                            double phi, double theta, double w_start, double w_end, double dw)
{
    std::vector<Vector3D> trajectory_s;
    std::vector<Vector3D> velocity_s;
    std::vector<double> spectrum_data;
    std::vector<double> spectrum_and_frequency;

    for (int i = 0; i < trajectory.size(); i+=3)
    {
        Vector3D pos(trajectory[i], trajectory[i + 1], trajectory[i + 2]);
        Vector3D vel(velocity[i], velocity[i + 1], velocity[i + 2]);

        trajectory_s.push_back(pos);
        velocity_s.push_back(vel);
    }
    std::span<const Vector3D> trajectory_span(trajectory_s);
    std::span<const Vector3D> velocity_span(velocity_s);
    std::span<const double> time_span(time);

    RadiationSpectrum spectrum(trajectory_span, velocity_span, time_span);
    spectrum_data = spectrum.calculateSpectrum(phi, theta, w_start, w_end, dw);

    for (int i = 0; i < spectrum_data.size(); i++)
    {
        spectrum_and_frequency.push_back(w_start + i*dw);
        spectrum_and_frequency.push_back(spectrum_data[i]);
    }

    return spectrum_and_frequency;
}


std::vector<std::vector<double>> color_map(std::vector<double> time, std::vector<double> trajectory, std::vector<double> velocity,
                                        double phi, double theta, double w_start, double w_end, double dw, double dtheta)
{
    
    double pi = std::numbers::pi;
    int steps_w = static_cast<int>((w_end - w_start) / dw) + 1;
    int steps_theta = static_cast<int>(2 * pi / dtheta) + 1;

    std::vector<Vector3D> trajectory_s;
    std::vector<Vector3D> velocity_s;
    std::vector<double> spectrum_data;
    std::vector<double> spectrum_and_frequency;
    std::vector<std::vector<double>> color_map(steps_theta, std::vector<double>(steps_w));

    for (int i = 0; i < trajectory.size(); i += 3)
    {
        Vector3D pos(trajectory[i], trajectory[i + 1], trajectory[i + 2]);
        Vector3D vel(velocity[i], velocity[i + 1], velocity[i + 2]);

        trajectory_s.push_back(pos);
        velocity_s.push_back(vel);
    }
    std::span<const Vector3D> trajectory_span(trajectory_s);
    std::span<const Vector3D> velocity_span(velocity_s);
    std::span<const double> time_span(time);

    RadiationSpectrum spectrum(trajectory_span, velocity_span, time_span);

    for (int i = 0; i < steps_theta; ++i)
    {
        double theta_map = i * 2 * pi / (steps_theta - 1);
        for (int j = 0; j < steps_w; ++j)
        {
            double w_map = j * w_end / (steps_w - 1);
            color_map[i][j] = spectrum.calculateIntensity(phi, theta_map, w_map);
        }
    }

    return color_map;
}