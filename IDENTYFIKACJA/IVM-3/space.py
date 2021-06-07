import random as r 
import numpy as np
from func import call

W = 0.5
c1 = 2.0  #1,5 lub 1,0
c2 = 2.0  # 1,5 lub 1,0

class Space():
    def __init__(self, target, epsilon, num_particles):
        self.target = target
        self.epsilon = epsilon
        self.num_particles = num_particles
        self.particles = []
        self.gbest_value = float('inf')
        self.gbest_position = np.array([
            r.uniform(0.0018, 0.0025),
            r.uniform(175.0, 177.0),
            r.uniform(19.4 * pow(10, 3), 19.6 * pow(10, 3)),
            r.uniform(0.0001 * 3 * pow(10, 10), 0.0002 * 3 * pow(10, 10)),
            r.uniform(150.00 * pow(10, 3), 152.00 * pow(10, 3)),
            r.uniform(0.97, 0.98),
            r.uniform(5.7, 5.8),
            r.uniform(0.8, 1.2),
            r.uniform(0.0, 0.1),
            0.0,
            r.uniform(0.0 * pow(10, 13), 0.1 * pow(10, 13)),
            r.uniform(0.0005 * pow(10, 13) ,0.0007 * pow(10, 13)),
            r.uniform(0.15, 0.2)
        ])
        
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
            p.velocity = new_velocity
            p.move()