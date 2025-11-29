class PID:
    def __init__(self, Kp, Ki, Kd, setpoint=0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.setpoint = setpoint
        self.integral = 0
        self.previous_error = 0

    def update(self, dt, measurement, rate=None):
        """Update the PID controller and return the control output.
        dt: time step in seconds
        measurement: observed process variable
        rate: optional rate of change of measurement
        """
        error = self.setpoint - measurement
        self.integral += error * dt
        derivative = 0.0
        if rate is not None:
            derivative = rate
        else:
            derivative = (error - self.previous_error) / dt if dt > 0 else 0

        output = (self.Kp * error) + (self.Ki * self.integral) + (self.Kd * derivative)

        self.previous_error = error

        return output
    
    def set_setpoint(self, setpoint):
        self.setpoint = setpoint
        self.integral = 0
        self.previous_error = 0

    def reset(self):
        self.integral = 0
        self.previous_error = 0

    def tune(self, Kp, Ki, Kd):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    pid = PID(Kp=5.0, Ki=0.3, Kd=0.05, setpoint=100)
    measurement = 0
    dt = 0.1
    outputs = [0]
    measurements = [measurement]
    clamp = 100.0
    for i in range(100):
        output = pid.update(dt, measurement)
        measurement += min(output, clamp) * dt + 5.0 * (np.random.random() - 0.5) - 2 # some noise and drift
        outputs.append(output)
        measurements.append(measurement)
    plt.figure()
    plt.plot(measurements, label='Measurement')
    plt.axhline(pid.setpoint, color='r', linestyle='--', label='Setpoint')
    plt.legend()
    plt.grid()
    plt.ylim(0, pid.setpoint + 40)

    plt.figure()
    plt.plot(outputs, label='PID Output')
    plt.axhline(clamp, color='r', linestyle='--', label='Clamp')
    plt.legend()
    plt.grid()
    plt.show()