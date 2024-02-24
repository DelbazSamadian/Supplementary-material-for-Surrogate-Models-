# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 09:46:22 2023

@author: User
"""

import numpy as np

def calculate_arias_intensity_and_cav(acceleration, time_step):
    """
    Calculate Arias Intensity and Cumulative Absolute Velocity (CAV) from a ground motion record.

    Parameters:
    - acceleration (numpy.ndarray): Ground motion acceleration time history.
    - time_step (float): Time step between acceleration data points.

    Returns:
    - arias_intensity (float): Arias Intensity value.
    - cav (float): Cumulative Absolute Velocity (CAV) value.
    """
    # Calculate Arias Intensity
    arias_intensity = np.trapz(np.square(acceleration), dx=time_step) * np.pi / 2.0

    # Calculate Cumulative Absolute Velocity (CAV)
    velocity = np.cumsum(np.abs(acceleration)) * time_step
    cav = np.max(velocity)

    return arias_intensity, cav

# Example usage
if __name__ == "__main__":
    # Sample ground motion acceleration data (replace with your own data)
    acceleration = np.array([0.1, 0.2, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3, -0.2])

    # Time step between acceleration data points (in seconds)
    time_step = 0.1

    # Calculate Arias Intensity and CAV
    arias_intensity, cav = calculate_arias_intensity_and_cav(acceleration, time_step)

    # Print the results
    print(f"Arias Intensity: {arias_intensity} m/s^2·s")
    print(f"CAV: {cav} m/s")

