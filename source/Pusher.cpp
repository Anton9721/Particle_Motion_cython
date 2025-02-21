#include "Pusher.hpp"

void PusherEuler::push(std::span<Particle> particles, const EMField &field, double t, double dt) const
{
    for (auto &particle : particles)
    {
        Vector3D e_field = field.get_electric_field(particle.position, t);
        Vector3D b_field = field.get_magnetic_field(particle.position, t);

        Vector3D force = (e_field + cross(particle.velocity, b_field)) * particle.charge;
        // force = PusherEuler::LorenzForce();

        Vector3D acceleration = force / particle.mass;

        particle.velocity += acceleration * dt;

        particle.position += particle.velocity * dt;
    }
}

Vector3D RungeKutta4::get_Lorentz_force(const Particle &particle, const EMField &field, const Vector3D &r, const Vector3D &velocity, double t) const
{
    Vector3D e_field = field.get_electric_field(r, t);
    Vector3D b_field = field.get_magnetic_field(r, t);

    return particle.charge * (e_field + cross(particle.velocity, b_field));
}

void RungeKutta4::push(std::span<Particle> particles, const EMField &field, double t, double dt) const
{
    for (auto &particle : particles)
    {

        Vector3D r0 = particle.position;
        Vector3D v0 = particle.velocity;
        double t0 = t;

        Vector3D k1_v = dt * get_Lorentz_force(particle, field, r0, v0, t0);
        Vector3D k1_r = dt * v0;

        Vector3D k2_v = dt * get_Lorentz_force(particle, field, r0 + 0.5 * k1_r, v0 + 0.5 * k1_v, t0 + 0.5 * dt);
        Vector3D k2_r = dt * (v0 + 0.5 * k1_v);

        Vector3D k3_v = dt * get_Lorentz_force(particle, field, r0 + 0.5 * k2_r, v0 + 0.5 * k2_v, t0 + 0.5 * dt);
        Vector3D k3_r = dt * (v0 + 0.5 * k2_v);

        Vector3D k4_v = dt * get_Lorentz_force(particle, field, r0 + k3_r, v0 + k3_v, t + dt);
        Vector3D k4_r = dt * (v0 + k3_v);

        particle.velocity = particle.velocity + (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v) / 6.0;
        particle.position = particle.position + (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r) / 6.0;
    }
}

void PusherBoris::push(std::span<Particle> particles, const EMField &field, double t, double dt) const
{
    for (auto &particle : particles)
    {
        Vector3D v = particle.velocity;
        Vector3D x_n = particle.position;
        double m = particle.mass;
        double q = particle.charge;

        double gamma = 1.0 / std::sqrt(1.0 - dot(v, v));
        Vector3D u_n = v * gamma;
        Vector3D x_n_half = x_n + u_n * dt / (2.0 * gamma);

        Vector3D e_field = field.get_electric_field(x_n_half, t);
        Vector3D b_field = field.get_magnetic_field(x_n_half, t);

        Vector3D u_minus = u_n + (q * dt * e_field) / (2.0 * m);
        Vector3D t_vec = b_field * (q * dt / (2.0 * m));
        Vector3D s = 2.0 * t_vec / (1.0 + dot(t_vec, t_vec));
        Vector3D u_plus = u_minus + cross(u_minus + cross(u_minus, t_vec), s);

        Vector3D u_n1 = u_plus + (q * dt * e_field) / (2.0 * m);
        double gamma_n1 = std::sqrt(1.0 + dot(u_n1, u_n1));
        Vector3D x_n1 = x_n_half + dt * u_n1 / (2.0 * gamma_n1);

        particle.velocity = u_n1 / gamma_n1;
        particle.position = x_n1;
    }
}

void PusherBorisRR::push(std::span<Particle> particles, const EMField &field, double t, double dt) const
{
    const double norm_const = -1.18 * 1e-8;

    for (auto &particle : particles)
    {
        Vector3D v = particle.velocity;
        Vector3D x_n = particle.position;
        double m = particle.mass;
        double q = particle.charge;

        double gamma = 1.0 / std::sqrt(1.0 - dot(v, v));
        Vector3D u_n = v * gamma;
        Vector3D x_n_half = x_n + u_n * dt / (2.0 * gamma);

        Vector3D e_field = field.get_electric_field(x_n_half, t);
        Vector3D b_field = field.get_magnetic_field(x_n_half, t);

        Vector3D lorentz_force = q * (e_field + cross(v, b_field));
        Vector3D rad_force = norm_const * gamma * (dot(lorentz_force, lorentz_force)
                                    - dot(v, lorentz_force) * dot(v, lorentz_force)) * v;

        Vector3D u_minus = u_n + (q * dt * e_field) / (2.0 * m);
        Vector3D t_vec = b_field * (q * dt / (2.0 * m));
        Vector3D s = 2.0 * t_vec / (1.0 + dot(t_vec, t_vec));
        Vector3D u_plus = u_minus + cross(u_minus + cross(u_minus, t_vec), s);

        Vector3D u_n1 = u_plus + (q * dt * e_field) / (2.0 * m) + rad_force * dt / m;
        double gamma_n1 = std::sqrt(1.0 + dot(u_n1, u_n1));
        Vector3D x_n1 = x_n_half + dt * u_n1 / (2.0 * gamma_n1);

        particle.velocity = u_n1 / gamma_n1;
        particle.position = x_n1;
    }
}



// void PusherBorisRR::push(std::span<Particle> particles, const EMField &field, double t, double dt) const
// {
//     const double norm_const = -1.18 * 1e-8;

//     for (auto &particle : particles)
//     {
//         Vector3D v = particle.velocity;
//         Vector3D x_n = particle.position;
//         double m = particle.mass;
//         double q = particle.charge;

//         double gamma = 1.0 / std::sqrt(1.0 - dot(v, v));
//         Vector3D u_n = v * gamma;
//         Vector3D x_n_half = x_n + u_n * dt / (2.0 * gamma);

//         Vector3D e_field = field.get_electric_field(x_n_half, t);
//         Vector3D b_field = field.get_magnetic_field(x_n_half, t);

//         Vector3D lorentz_force = q * (e_field + cross(v, b_field));

//         // Vector3D rad_force = norm_const *
//         //                  (cross(lorentz_force, b_field) - (dot(v, e_field))*e_field +
//         //                   gamma * gamma * (dot(lorentz_force, lorentz_force) - (dot(v, e_field)) * (dot(v, e_field))) * v);
//         Vector3D rad_force = norm_const * gamma * (dot(lorentz_force, lorentz_force) 
//                             - dot(v, lorentz_force) * dot(v, lorentz_force)) * v;

//         Vector3D u_minus = u_n + dt * (q * e_field + rad_force) / (2.0 * m);

//         Vector3D t_vec = b_field * (q * dt / (2.0 * m * gamma));
//         Vector3D s = 2.0 * t_vec / (1.0 + dot(t_vec, t_vec));
//         Vector3D u_plus = u_minus + cross(u_minus + cross(u_minus, t_vec), s);

//         Vector3D u_n1 = u_plus + dt * (q * e_field + rad_force) / (2.0 * m);
//         double gamma_n1 = std::sqrt(1.0 + dot(u_n1, u_n1));
//         Vector3D x_n1 = x_n_half + dt * u_n1 / (2.0 * gamma_n1);

//         particle.velocity = u_n1 / gamma_n1;
//         particle.position = x_n1;
//     }
// }


