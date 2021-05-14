import random as r 
import numpy as np
from func import square_error as call

W = 0.5
c1 = 2.1
c2 = 5.1

class Space():
    def __init__(self, target, epsilon, num_particles):
        self.target = target
        self.epsilon = epsilon
        self.num_particles = num_particles
        self.particles = []
        self.gbest_value = np.longdouble('inf')
        self.gbest_position = np.array([r.uniform(1,1000), r.uniform(0,1), r.uniform(0,1), r.uniform(1, 10_000), r.uniform(0,1), r.uniform(1, 90_000), r.uniform(0,1)], dtype=np.longdouble)
        
    def print_particles(self):
        for p in self.particles:
            p.__str__()
    
    def fitness(self, particle):
        return call(particle.position)
    
    def set_pbest(self):
        for p in self.particles:
            fitness_candidate = self.fitness(p)
            if(p.pbest_value > fitness_candidate):
                p.pbest_value = fitness_candidate
                p.pbest_position = p.position
    
    def set_gbest(self):
        for p in self.particles:
            best_fitness_candidate = self.fitness(p)
            if(self.gbest_value > best_fitness_candidate):
                self.gbest_value = best_fitness_candidate
                self.gbest_position = p.position
              
    def move_particles(self):
        W = r.random()
                
        for p in self.particles:
            new_velocity = W*p.velocity + c1*r.random()*(p.pbest_position-p.position) + r.random()*c2*(self.gbest_position - p.position)
            
            for i in range(len(new_velocity)):
                if new_velocity[i] > dict_constraints[i][1] or new_velocity[i] < dict_constraints[i][0]:
                    new_velocity[i] *= 0.95
            p.velocity = new_velocity
            p.move()
        
        

dict_constraints = {
    0: [1.0, 1000.0],
    1: [0.0, 1.0],
    2: [0.0, 1.0],
    3: [1.0, 10_000.0],
    4: [0.0, 1.0],
    5: [1.0, 90_000.0],
    6: [0.0, 1.0]
}            