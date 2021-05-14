import random as r
import numpy as np

class Particle():
    def __init__(self):
        self.position = np.array([r.uniform(5.5,6.5), r.uniform(6.5,8.5)])
        self.pbest_position = self.position
        self.pbest_value = float('inf')
        self.velocity = np.zeros((2))
        
    def __str__(self):
        print("Position ", self.position, " pbest: ", self.pbest_position)
    
    def move(self):
        self.position = self.position + self.velocity