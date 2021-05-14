import random as r
import numpy as np

class Particle():
    def __init__(self):
        self.position = np.array([r.uniform(490, 510), r.uniform(0.4, 0.6), r.uniform(0.4, 0.6), r.uniform(4900, 5100), r.uniform(0.4, 0.6), r.uniform(44000, 46000), r.uniform(0.4, 0.6)], dtype=np.longdouble)
        self.pbest_position = self.position
        self.pbest_value = np.longdouble('inf')
        self.velocity = np.zeros((7), dtype=np.longdouble)
        
    def __str__(self):
        print("Position ", self.position, " pbest: ", self.pbest_position)
    
    def move(self):
        self.position = self.position + self.velocity