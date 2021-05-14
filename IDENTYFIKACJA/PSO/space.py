import random as r 
import numpy as np
from func import call

W = 0.5
c1 = 2.1
c2 = 5.1

class Space():
    def __init__(self, target, epsilon, num_particles, func):
        self.target = target
        self.epsilon = epsilon
        self.num_particles = num_particles
        self.particles = []
        self.gbest_value = float('inf')
        self.gbest_position = np.array([r.uniform(5.75,6.95), r.uniform(6.5,8.5)])
        self.func = func
        
    def print_particles(self):
        for p in self.particles:
            p.__str__()
    
    def fitness(self, particle):
        return call(particle.position, self.func)
    
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
            