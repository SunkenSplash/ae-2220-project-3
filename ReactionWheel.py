import numpy as np

class ReactionWheel:
    def __init__(self, max_torque: float, inertia: float, axis: np.array):
        self.min_torque = -max_torque
        self.max_torque = max_torque
        self.inertia = inertia
        self.axis = axis
        self.angular_velocity = 0.0
        # Integrated wheel angle (radians). This is the wheel's rotation
        # angle in the wheel frame and is updated in `update()` so callers
        # can read a continuously integrated angle without re-integrating
        # histories themselves.
        self.angle = 0.0
        self.torque = 0.0

    def set_torque(self, torque: float):
        self.torque = np.clip(torque, self.min_torque, self.max_torque)

    def update(self, dt: float):
        # angular acceleration = torque / inertia
        # negative because of equal-and-opposite reaction
        self.angular_velocity -= (self.torque / self.inertia) * dt
        # Integrate the wheel angle using the (new) angular velocity so the
        # ReactionWheel maintains its own angle state.
        self.angle += self.angular_velocity * dt
        self.angle = self.angle