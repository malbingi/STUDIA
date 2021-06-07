import particle as p
import space as s
from saveToFile import saveToFile
from calculations_for_plots import ivm_calculations as calculations

optimization_target = 0.0
epsilon = 0.1


def main():
    search_space = s.Space(optimization_target, epsilon, 100)
    particles_vector = [p.Particle()
                        for _ in range(search_space.num_particles)]
    search_space.particles = particles_vector
    search_space.set_pbest()
    search_space.set_gbest()

    iteration = 0
    n_iterations = 1000
    val_prev = val_prev_1 = 0.0
    while iteration < n_iterations or search_space.gbest_value < epsilon:
        print(iteration)
        search_space.set_pbest()
        search_space.set_gbest()

        if iteration == 0:
            val_prev = search_space.gbest_value
        else:
            val_prev_1 = val_prev
            val_prev = search_space.gbest_value

        #print("POŁOŻENIE: ", search_space.gbest_position, " WARTOŚĆ: ", search_space.gbest_value)

        saveToFile(search_space.gbest_value, "main_18_100.csv")
        saveToFile(search_space.gbest_position, "coordinates_18_100.csv")

        if iteration > 10:
            if abs(search_space.gbest_value - search_space.target) <= search_space.epsilon or val_prev_1 == val_prev:
                break

        search_space.move_particles()
        iteration += 1
    #print("Solution: ", search_space.gbest_position, " value ", search_space.gbest_value)

def calculate_based_on_vector():
    vector = [2.45212129e-03,  1.76819149e+02,  1.93909885e+04,  1.20589593e+07,
              1.38080807e+05,  9.73000000e-01,  5.77000000e+00,  6.36456556e-01,
              - 5.70121100e-02, - 5.13174918e-01,  0.00000000e+00, 8.82004319e+09,
              9.90984827e-03]
    dict_T_epsilon_dot = {
        0: [950.0, 0.1, 10.0],
        1: [1050.0, 0.1, 10.0],
        2: [1200.0, 0.1, 10.0],
        3: [900.0, 1.0, 1.0],
        4: [1050.0, 1.0, 1.0],
        5: [1200.0, 1.0, 1.0],
        6: [900.0, 10.0, 0.1],
        7: [1050.0, 10.0, 0.1],
        8: [1200.0, 10.0, 0.1]
    }
    dict_names = {
        0: '950_01.csv',
        1: '1050_01.csv',
        2: '1200_01.csv',
        3: '900_1.csv',
        4: '1050_1.csv',
        5: '1200_1.csv',
        6: '900_10.csv',
        7: '1050_10.csv',
        8: '1200_10.csv',
    }
    for i in range(9):
        temp, epsilon_dot, time = dict_T_epsilon_dot[i]
        name = dict_names[i]
        calculations(vector, temp, epsilon_dot, time, name)
    


if __name__ == "__main__":
    #main()
    calculate_based_on_vector()
