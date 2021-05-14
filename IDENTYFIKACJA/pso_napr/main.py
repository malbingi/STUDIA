import particle as p
import space as s
import queue as q
from saveToFile import saveToFile as save
from func import sigma_p as sigma

epsilon = 0.01

dict_T_epsilon_dot = {
    0: [850, 0.1],
    1: [850, 1],
    2: [850, 10],
    3: [1000, 0.1],
    4: [1000, 1],
    5: [1000, 10],
    6: [1150, 0.1],
    7: [1150, 1],
    8: [1150, 10]
}

def main():
    search_space = s.Space(0.0, epsilon, 500)
    particles_vector = [p.Particle() for _ in range(search_space.num_particles)]
    search_space.particles = particles_vector
    search_space.set_pbest()   
    search_space.set_gbest()     
    
    iteration = 0
    n_iterations = 100
    val_prev = val_prev_1 = 0.0
    while iteration < n_iterations or search_space.gbest_value <= epsilon:        
        search_space.set_pbest()
        search_space.set_gbest()
        
        # JEŚLI PRZEZ TRZY ITERACJE NIE ZMIENI SIE WARTOŚĆ => BREAK
         
        if iteration == 0:
            val_prev = search_space.gbest_value
        else: 
            val_prev_1 = val_prev
            val_prev = search_space.gbest_value
        
        #print("POŁOŻENIE: ", search_space.gbest_position, " WARTOŚĆ: ", search_space.gbest_value)
        if abs(search_space.gbest_value - search_space.target) <= search_space.epsilon or val_prev_1 == val_prev:
            break
        
        search_space.move_particles()
        iteration += 1
        save(search_space.gbest_value, 0)
    print("Solution: ", search_space.gbest_position, " value ", search_space.gbest_value)
    file = open("model.txt", 'a')
    file.write(str(search_space.gbest_position) + " \n")
    file.close()
    return search_space.gbest_position
    
if __name__ == "__main__":
    vector = main()
    for i in range(9):
        T, epsilon_dot = dict_T_epsilon_dot[i]
        for j in range(21):
           save(sigma(vector, 0.05*(j+1), epsilon_dot, T), i+1)