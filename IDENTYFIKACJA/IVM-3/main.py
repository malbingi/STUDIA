import particle as p
import space as s
from saveToFile import saveToFile

optimization_target = 0.0
epsilon = 0.1

def main():
    search_space = s.Space(optimization_target, epsilon, 20)
    particles_vector = [p.Particle() for _ in range(search_space.num_particles)]
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
        
        saveToFile(search_space.gbest_value, "main12.csv")
        saveToFile(search_space.gbest_position, "coordinates12.csv")
        
        if iteration > 10:
            if abs(search_space.gbest_value - search_space.target) <= search_space.epsilon or val_prev_1 == val_prev:
                break
        
        search_space.move_particles()
        iteration += 1
    #print("Solution: ", search_space.gbest_position, " value ", search_space.gbest_value)
    
if __name__ == "__main__":
    main()
