#include "Writer.hpp"

void WriterConsole::write(std::span<Particle> particles, double t)
{
    for (const auto &particle : particles)
    {
        std::cout << "Particle position: ("
                  << particle.position.x << ", "
                  << particle.position.y << ", "
                  << particle.position.z << ") \t";

        std::cout << "Particle velocity: ("
                  << particle.velocity.x << ", "
                  << particle.velocity.y << ", "
                  << particle.velocity.z << ")\n\n";
    }
}

WriterTXT::WriterTXT(std::string &filename)
{
    file = std::ofstream(filename);
}

void WriterTXT::write(std::span<Particle> particles, double t)
{

    file << t << std::setw(2);
    for (const auto &particle : particles)
    {
        file << particle.position << std::setw(2) << particle.velocity;
    }
    file << "\n";
}

std::vector<double> WriterPyhton::writer_numpy(std::span<Particle> particles, double t)
{
    std::vector<double> data_result_part;
    for (const auto &particle : particles)
    {
        data_result_part.push_back(t);
        data_result_part.push_back(particle.position.x);
        data_result_part.push_back(particle.position.y);
        data_result_part.push_back(particle.position.z);
        data_result_part.push_back(particle.velocity.x);
        data_result_part.push_back(particle.velocity.y);
        data_result_part.push_back(particle.velocity.z);
    }


    return data_result_part;
}

WriterSpectrum::WriterSpectrum(std::string &filename)
{
    file = std::ofstream(filename);
}

void WriterSpectrum::write(std::vector<double> spectrum_data, double w_start, double w_end, double dw)
{
    double w = w_start;
    for (const auto &row : spectrum_data)
    {
        file << w << " " << row << "\n";
        w += dw;
    }
}

TrajectoryReader::TrajectoryReader(const std::string &filename)
{
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double time;
        char dummy;
        Vector3D position, velocity;

        if (iss >> time >> dummy >> position.x >> dummy >> position.y >> dummy >> position.z >> dummy >>
            dummy >> velocity.x >> dummy >> velocity.y >> dummy >> velocity.z >> dummy)
        {
            trajectory_.push_back(position);
            velocity_.push_back(velocity);
            time_.push_back(time);
        }
    }

    file.close();
}




WriterSpectrumMap::WriterSpectrumMap(std::string& filename)
{
    file = std::ofstream(filename);
}

void WriterSpectrumMap::write(std::vector<std::vector<double>> colour_map)
{
    for (const auto &row : colour_map)
    {
        for (double value : row)
        {
            file << value << " ";
        }
        file << "\n";
    }
    file.close();
}
