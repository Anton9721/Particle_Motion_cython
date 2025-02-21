#pragma once
#include <span>
#include "Particle.hpp"
#include "RadiationSpectrum.hpp"
#include <string>
#include <span>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

class Writer
{
public:
    virtual void write(std::span<Particle> particles, double t) = 0;
    virtual ~Writer() = default;
};

// Вывод в консоль
class WriterConsole : public Writer
{
public:
    void write(std::span<Particle> particles, double t) override;
};

// Вывод в numpy массив
class WriterPyhton
{
public:
    virtual ~WriterPyhton() {}

    std::vector<double> writer_numpy(std::span<Particle> particles, double t);
};

// Запись в файл txt времени, позиции и скорости
class WriterTXT : public Writer
{
private:
    std::ofstream file;

public:
    WriterTXT(std::string &filename);
    virtual ~WriterTXT() {}
    void write(std::span<Particle> particles, double t) override;
};

// Запись в файл txt спектра
class WriterSpectrum
{
private:
    std::ofstream file;

public:
    WriterSpectrum(std::string &filename);
    virtual ~WriterSpectrum() {}
    void write(std::vector<double> data_spectrum, double w_start, double w_end, double dw);
};

// Запись в файл txt спектра для цветовой карты
class WriterSpectrumMap
{
private:
    std::ofstream file;

public:
    WriterSpectrumMap(std::string& filename);
    virtual ~WriterSpectrumMap() {}
    void write(std::vector<std::vector<double>> colour_map);
};

// Считывание данных из файла txt с траекториями и скоростями
class TrajectoryReader
{
public:
    TrajectoryReader(const std::string &filename);

    std::span<const Vector3D> getTrajectory() const { return trajectory_; }
    std::span<const Vector3D> getVelocity() const { return velocity_; }
    std::span<const double> getTime() const { return time_; }

private:
    std::vector<Vector3D> trajectory_;
    std::vector<Vector3D> velocity_;
    std::vector<double> time_;
};
