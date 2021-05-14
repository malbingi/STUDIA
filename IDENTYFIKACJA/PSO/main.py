import particle as p
import space as s
import queue as q
from saveToFile import saveToFile

epsilon = 0.0000001
# WYBÓR FUNKCJI
__function__ = 2

def main():
    search_space = s.Space(0.0, epsilon, 15, __function__)
    particles_vector = [p.Particle() for _ in range(search_space.num_particles)]
    search_space.particles = particles_vector
    search_space.set_pbest()   
    search_space.set_gbest()     
    
    iteration = 0
    n_iterations = 1000
    val_prev = val_prev_1 = 0.0
    while iteration < n_iterations or search_space.gbest_value < epsilon:        
        search_space.set_pbest()
        search_space.set_gbest()
        
        # JEŚLI PRZEZ TRZY ITERACJE NIE ZMIENI SIE WARTOŚĆ => BREAK
        if iteration == 0:
            val_prev = search_space.gbest_value
        else: 
            val_prev_1 = val_prev
            val_prev = search_space.gbest_value
        
        print("POŁOŻENIE: ", search_space.gbest_position, " WARTOŚĆ: ", search_space.gbest_value)
        saveToFile(__function__, search_space.gbest_value)
        if abs(search_space.gbest_value - search_space.target) <= search_space.epsilon or val_prev_1 == val_prev:
            break
        
        search_space.move_particles()
        iteration += 1
    print("Solution: ", search_space.gbest_position, " value ", search_space.gbest_value)
    
if __name__ == "__main__":
    main()