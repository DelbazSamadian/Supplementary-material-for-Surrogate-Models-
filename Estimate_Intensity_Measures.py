# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:37:28 2023

@author: User
"""

import numpy as np
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


def duration_params(acceleration):
    """
    Calculate D95-D5 from a ground motion record.

    Parameters:
    - acceleration (numpy.ndarray): Ground motion acceleration time history.

    Returns:
    - d95_d5 (float): D95-D5 value.
    """
    abs_acceleration = np.abs(acceleration)
    d95 = np.percentile(abs_acceleration, 95)
    d5 = np.percentile(abs_acceleration, 5)
    d95_d5 = d95 - d5
    return d95_d5

def main():
    # Load the ground motion record
    file_path = "Time_History.txt"
    try:
        data = np.loadtxt(file_path)
    except Exception as e:
        print(f"Error loading ground motion record from {file_path}: {e}")
        return

    time_values = data[:, 0]
    acceleration_values = data[:, 1]

    # Calculate CAV, Ia, and D95-D5 for the single record
    arias_intensity, cav = calculate_arias_intensity_and_cav (acceleration_values, np.mean(np.diff(time_values)))
    d95_d5 = duration_params(acceleration_values)

    print("Arias Intensity:", arias_intensity)
    print("CAV:", cav)
    print("D95-D5:", d95_d5)

if __name__ == "__main__":
    main()
