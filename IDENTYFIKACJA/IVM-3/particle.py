import random as r
import numpy as np


class Particle():
    def __init__(self):
        self.position = np.array([
            r.uniform(0.002, 0.0025),
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
            r.uniform(0.1, 0.2)
        ])
        self.pbest_position = self.position
        self.pbest_value = float('inf')
        self.velocity = np.zeros((13))

    def __str__(self):
        print("Position ", self.position, " pbest: ", self.pbest_position)

    def move(self):
        self.position = self.position + self.velocity
        self.limit_vector()

    def limit_vector(self):
        if self.position[0] < 0.0:
            self.position[0] *= -1.
        if self.position[0] > 0.01:
            self.position[0] *= 0.5

        if self.position[1] < 0.0:
            self.position[1] *= -1.

        if self.position[1] < 10.:
            self.position[1] += 10.
        if self.position[1] > 500.:
            self.position[1] *= 0.5

        if self.position[2] < 0.0:
            self.position[2] *= -1.
        if self.position[2] < 1000.0:
            self.position[2] += 1000.
        if self.position[2] > 100_000.0:
            self.position[2] *= 0.5

        if self.position[3] < 0.0:
            self.position[3] *= -1.
        if self.position[3] > 3000000000.0:
            self.position[3] *= 0.5

        if self.position[4] < 0.0:
            self.position[4] *= -1.
        if self.position[4] < 10_000.0:
            self.position[4] += 10_000.0
        if self.position[4] > 500_000.0:
            self.position[4] *= 0.5

        if self.position[5] != 0.973:
            self.position[5] = 0.973

        if self.position[6] != 5.77:
            self.position[6] = 5.77

        if self.position[7] < 0.0:
            self.position[7] *= -1.
        if self.position[7] < 0.9:
            self.position[7] *= 1.5
        if self.position[7] > 1.0:
            self.position[7] *= 0.25

        if self.position[8] < 0.0:
            self.position[8] *= -1.
        if self.position[8] > 10.0:
            self.position[8] *= 0.5

        if self.position[9] < 0.0:
            self.position[9] *= -1.
        if self.position[9] > 1.0:
            self.position[9] *= 0.5

        if self.position[10] != 0:
            self.position[10] = 0.0

        if self.position[11] < 0.0:
            self.position[11] *= -1.
        if self.position[11] > 1_000_000_000_000.0:
            self.position[11] -= 1_000_000_000_000.0

        if self.position[12] < 0.0:
            self.position[12] *= -1.
        if self.position[12] > 1.0:
            self.position[12] *= 0.5
