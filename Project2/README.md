# Project 2 — Lateral Vehicle Model

## Overview
This project builds a lateral dynamic model of a race car in Simulink and MATLAB. The workflow: `init.m` defines vehicle parameters/initial conditions, `gentrack.m` generates a closed-course centerline and borders, the Simulink model (`Project_2_Kinematic_Model.slx`) runs the vehicle with a driver, lateral dynamics, and coordinate transformation, and `animate.m` animates the lap and reports basic stats.

## Files
- `run.m` — entry point; calls `init.m`, `gentrack.m`, then launches the Simulink model and animation.
- `init.m` — sets vehicle mass, inertia, tire stiffness, geometry, desired speed, and initial states; creates `carDataBus` for Simulink.
- `gentrack.m` — creates `path` struct with center, inner, and outer waypoints for a 900 m straight × 200 m radius oval (15 m width).
- `Project_2_Kinematic_Model.slx` — Simulink model (driver, lateral dynamics, coordinate transform, XY graph, To Workspace blocks) to be completed/tuned.
- `animate.m` — runs the Simulink sim, animates the car polygon over the track, and computes simple race statistics.
- `raceStat.m` — post-processing utility used by `animate.m` to compute lap counts, speed, and off-track events.
- `Project2_Description_and_deliverbles.pdf` — assignment prompt with week-by-week requirements.

## How to Run
1) Open MATLAB/Simulink (R2023b+ recommended). Ensure the current folder is this `Project2` directory.
2) Run `run.m` in MATLAB. It will generate the track, execute the Simulink model, and animate results. If Simulink signals differ from `simout.X/Y/psi`, update `animate.m` or the `To Workspace` block names accordingly.
3) Adjust vehicle/track parameters in `init.m` and the target-speed/driver gains inside the Simulink model to hit lap-time goals without leaving the 15 m width.

## Notes / Pending Work
- Driver, Lateral Dynamics, and Transform subsystems in the Simulink model need completion and tuning (see `TODO.md`).
- Add/verify `To Workspace` exports for X, Y, psi, velocity, and time to feed `animate.m` and `raceStat`.
