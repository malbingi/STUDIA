import particle as p
import space as s
from saveToFile import saveToFile
from calculations_for_plots import ivm_calculations as calculations

optimization_target = 0.0
epsilon = 0.1


def main():
    search_space = s.Space(optimization_target, epsilon, 20)
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

        saveToFile(search_space.gbest_value, "main_19_20.csv")
        saveToFile(search_space.gbest_position, "coordinates_19_20.csv")

        if iteration > 10:
            if abs(search_space.gbest_value - search_space.target) <= search_space.epsilon or val_prev_1 == val_prev:
                break

        search_space.move_particles()
        iteration += 1
    #print("Solution: ", search_space.gbest_position, " value ", search_space.gbest_value)
    return search_space.gbest_position

def calculate_based_on_vector(vector):
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
    vector = main()
    calculate_based_on_vector(vector)
