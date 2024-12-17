import numpy as np
cimport numpy as cnp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "include/Library.hpp":
    vector[double] find_trajectory(
        string solver_name, string field_name,
        vector[double] position, vector[double] velocity,
        double mass, double charge,
        vector[double] electric_field, vector[double] magnetic_field,
        double time_start, double time_end, double dt, double omega_f, double e_initial_phase, double b_initial_phase
    )

    vector[double] spectrum(vector[double] time, vector[double] trajectory, 
        vector[double] velocity,  double phi, double theta, double w_start,
        double w_end, double dw)

    vector[vector[double]] color_map(vector[double] time, vector[double] trajectory, vector[double] velocity,
                                           double phi, double theta, double w_start, double w_end, double dw, double dtheta);


cdef vector[double] python_list_to_vector(cnp.ndarray[cnp.double_t, ndim=1] np_array):
    cdef vector[double] cpp_vector
    for i in range(np_array.shape[0]):
        cpp_vector.push_back(np_array[i])
    return cpp_vector

def vector_to_python_list(vector[double] cpp_vector):
    return [cpp_vector[i] for i in range(cpp_vector.size())]

def vector_vector_to_python_list(vector[vector[double]] cpp_vector):
    return [cpp_vector[i] for i in range(cpp_vector.size())]

def convert_to_numpy_array(vector[vector[double]] cpp_vector):
    rows = len(cpp_vector)
    cols = len(cpp_vector[0]) if rows > 0 else 0
    numpy_array = np.zeros((rows, cols), dtype=np.float64)
    for i in range(rows):
        for j in range(cols):
            numpy_array[i, j] = cpp_vector[i][j]
    
    return numpy_array

def trajectory_array(solver_name: str, field_name: str,
                     position, velocity,
                     mass: float, charge: float,
                     electric_field, magnetic_field,
                     time_start: float, time_end: float, dt: float, omega_f = 1.0,
                     e_initial_phase = 0.0, b_initial_phase = 0.0):


    all_solvers = ['PusherEuler', 'RungeKutta4', 'PusherBoris', 'PusherBorisRR']
    all_field = ['CrossEMField', 'GaussEMField', 'OscillatingEMField']

    if (solver_name not in all_solvers):
        print(f"Не существует метода решения {solver_name}. Проверьте правильность написания !")
        print(f"Допустимые решатели: {all_solvers}")
        return 0

    if field_name not in all_field:
        print(f"Не существует такого типа поля {solver_name}. Проверьте правильность написания !")
        print(f"Допустимые поля: {all_field}")
        return 0
    

   
    if isinstance(position, list):
        position = np.array(position, dtype=np.float64)
    if isinstance(velocity, list):
        velocity = np.array(velocity, dtype=np.float64)
    if isinstance(electric_field, list):
        electric_field = np.array(electric_field, dtype=np.float64)
    if isinstance(magnetic_field, list):
        magnetic_field = np.array(magnetic_field, dtype=np.float64)

    cdef string cpp_solver_name = solver_name.encode('utf-8')
    cdef string cpp_field_name = field_name.encode('utf-8')

    cdef vector[double] cpp_position = python_list_to_vector(position)
    cdef vector[double] cpp_velocity = python_list_to_vector(velocity)
    cdef vector[double] cpp_electric_field = python_list_to_vector(electric_field)
    cdef vector[double] cpp_magnetic_field = python_list_to_vector(magnetic_field)

    cdef vector[double] result_vector = find_trajectory(
        cpp_solver_name,
        cpp_field_name,
        cpp_position,
        cpp_velocity,
        mass,
        charge,
        cpp_electric_field,
        cpp_magnetic_field,
        time_start,
        time_end,
        dt, 
        omega_f, 
        e_initial_phase,
        b_initial_phase
    )

    return vector_to_python_list(result_vector)

def spectrum_array(data_trajectory,  phi: float, 
    theta: float, w_start: float, w_end: float, dw: float):
    
    time, trajectory, velocity = parse_data_trajectory(data_trajectory)
    
    cdef vector[double] cpp_time = python_list_to_vector(time)
    cdef vector[double] cpp_trajectory = python_list_to_vector(trajectory)
    cdef vector[double] cpp_velocity = python_list_to_vector(velocity)

    cdef vector[double] result_vector = spectrum(
        cpp_time, cpp_trajectory, cpp_velocity, phi, theta, w_start,
        w_end, dw
    )

    return vector_to_python_list(result_vector)


def color_map_array(data_trajectory,  phi: float, 
    theta: float, w_start: float, w_end: float, dw: float, dtheta: float):
    
    time, trajectory, velocity = parse_data_trajectory(data_trajectory)
    
    cdef vector[double] cpp_time = python_list_to_vector(time)
    cdef vector[double] cpp_trajectory = python_list_to_vector(trajectory)
    cdef vector[double] cpp_velocity = python_list_to_vector(velocity)

    cdef vector[vector[double]] result_vector = color_map(
        cpp_time, cpp_trajectory, cpp_velocity, phi, theta, w_start,
        w_end, dw, dtheta
    )

    return convert_to_numpy_array(result_vector)



# Рисование графиков
def parse_data_trajectory(data):
    time = []
    position = []
    velocities = []

    for i in range(0, len(data), 7): 
        time.append(data[i])  
        position.append(data[i+1])
        position.append(data[i+2])  
        position.append(data[i+3])  
        velocities.append(data[i+4])  
        velocities.append(data[i+5])
        velocities.append(data[i+6])

    return np.array(time), np.array(position), np.array(velocities)

def parse_data_spectrum(data):
    frequency = []
    intensity = []

    for i in range(0, len(data), 2): 
        frequency.append(data[i])  
        intensity.append(data[i+1])

    return np.array(frequency), np.array(intensity)

def kinetic_energy(velocities, mass):
    velocities = np.array(velocities).reshape(-1, 3)
    speeds_squared = np.sum(velocities**2, axis=1)
    return 0.5 * mass * speeds_squared


def writer_trajectory(data_trajectory, mass = 1):
    times, positions, velocities = parse_data_trajectory(data_trajectory)

    kinetic_energies = kinetic_energy(velocities, mass)

    fig = plt.figure(figsize=(12, 10))

    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot(positions[0::3], positions[1::3], positions[2::3], label="Trajectory")  # разбиение массива на 3 компонента
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title("3D Trajectory")
    ax1.legend()

    # График XY плоскости
    ax2 = fig.add_subplot(222)
    ax2.plot(positions[0::3], positions[1::3], color='g')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_title("XY Plane")

    # График XZ плоскости
    ax3 = fig.add_subplot(223)
    ax3.plot(positions[0::3], positions[2::3], color='r')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Z')
    ax3.set_title("XZ Plane")
    ax3.set_aspect('equal')

    # График YZ плоскости
    ax4 = fig.add_subplot(224)
    ax4.plot(positions[1::3], positions[2::3], color='b')
    ax4.set_xlabel('Y')
    ax4.set_ylabel('Z')
    ax4.set_title("YZ Plane")

    plt.figure()
    plt.plot(times, kinetic_energies, label="Kinetic Energy")
    plt.xlabel("Time")
    plt.ylabel("Kinetic Energy")
    plt.title("Kinetic Energy Over Time")
    plt.legend()

    plt.show()

def writer_spectrum(spectrum):
    frequencies, spectrum_values = parse_data_spectrum(spectrum)
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, spectrum_values, linestyle='-', color='b')

    plt.title("Зависимость спектра от частоты")
    plt.grid(True)

    plt.show() 

def writer_color_map(color_map, w_max):
    num_w = color_map.shape[1]
    num_theta = color_map.shape[0]
    w = np.linspace(0.1, w_max, num_w)  
    theta = np.linspace(0, 2 * np.pi, num_theta) 

    W, Theta = np.meshgrid(w, theta)

    color_map_log = np.log10(color_map)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
    contour = ax.contourf(Theta, W, color_map_log, levels=40, cmap="viridis") 

    cbar = fig.colorbar(contour, ax=ax, label="Logarithmic Intensity")
    cbar.set_label("Intensity (log scale)")

    plt.show()

