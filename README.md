# Quadrotor GNSS-Style Navigation Using a MOCAP Pseudolite System

This repository contains MATLAB code and experimental datasets used to estimate a quadrotor drone’s 3D trajectory using **pseudorange measurements** generated from a motion-capture (MOCAP) system configured to emulate a GNSS-like pseudolite network.

A detailed summary of the methods, theory, and results can be found in the project poster:

📄 **[Download the Project Poster (PDF)](poster.pdf)**  

---

## 📂 Repository Structure
.
├── Calibration_0_7B_Space.mat      # Static pseudorange measurements 
├── Calibration_1_7B_Space.mat      # Static pseudorange measurements
├── Calibration_2_7B_Space.mat      # Static pseudorange measurements
├── Calibration_3_7B_Space.mat      # Static pseudorange measurements
├── Calibration_4_7B_Space.mat      # Static pseudorange measurements
├── Error_Ellipse_Data.mat          # Additional uncertainty data
├── flight_data4_7B_Space.mat       # Main recorded MOCAP pseudorange and satellite location data for test flight
├── GPS_Nav_3D.m                    # Main MATLAB script: processing, estimation, plots
└── README.md                       # Documentation
