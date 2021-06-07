import random as r
import numpy as np


class Particle():
    def __init__(self):
        self.position = np.array([
            r.uniform(0.0018, 0.0025),
            r.uniform(175.0, 177.0),
            r.uniform(19.4 * pow(10, 3), 19.6 * pow(10, 3)),
            r.uniform(0.0001 * 3 * pow(10, 10), 0.0002 * 3 * pow(10, 10)),
            r.uniform(150.00 * pow(10, 3), 152.00 * pow(10, 3)),
            0.973,
            5.77,
            r.uniform(0.8, 1.2),
            r.uniform(0.0, 0.1),
            r.uniform(0.2, 0.3),
            0,
            r.uniform(0.0005 * pow(10, 13), 0.0007 * pow(10, 13)),
            r.uniform(0.15, 0.2)
        ])
        self.pbest_position = self.position
        self.pbest_value = float('inf')
        self.velocity = np.zeros((13))

    def __str__(self):
        print("Position ", self.position, " pbest: ", self.pbest_position)

    def move(self):
        self.position = self.position + self.velocity