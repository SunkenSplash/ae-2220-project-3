import numpy as np

class Thruster:
    def __init__(self, max_thrust: float, position: np.array, angle: float):
        self.max_thrust = max_thrust
        self.thrust = 0.0
        self.position = position
        self.angle = angle

    def set_thrust(self, thrust: float):
        self.thrust = np.clip(thrust, 0.0, self.max_thrust)