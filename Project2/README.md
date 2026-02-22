# Project 2 — Week 2

## Team 4
Dalys Guajardo, Juan Lopez, Ian Wilhite

## Summary

This project builds a lateral dynamic vehicle model in Simulink to drive around an oval race track following a set speed of 15 m/s without leaving the track boundaries.

**Week 1** setup the track geometry and animation. 

**Week 2** completes lateral dynamics model & ran closed loop simulation. The Simulink model contains four subsystems connected in a feedback loop:

1. **Driver Model** — a MATLAB Function (`SimplestDriver`) that takes vehicle position (X, Y), heading (psi), and yaw rate (omega) and computes the front steering angle delta_f using a PD controller with Ackermann feedforward (FWD). The reference heading comes from track geometry. Gains Kp = 0.5 and Kd = 2.0 were selected to provide an overdamped smooth motion.

![](./figures/driver_model.png)

2. **Lateral Dynamics** — tire slip angles feed a linear cornering-stiffness model (C_α = 40,000 N/rad, F_max = 2,793 N) to produce lateral tire forces, which drive kinematic equations for lateral velocity (v_y) and yaw rate (omega).

![](./figures/kinematics_model.png)

3. **Transformation/Rotation** — integrates the body-frame velocities rotated by heading angle psi to produce world-frame X, Y velocities & coordinates.

4. **XY Graph** — plots vehicle path.

X, Y, and psi are fed back to the Driver Model. The target speed was set to **15 m/s**. Maximum lateral acceleration with these tire parameters is 5.6 m/s², giving a curve speed limit of √(5.6 × 200) ≈ 33 m/s; 15 m/s provides a large margin and the vehicle stays within track boundaries for the full lap. Notably, it hugs the outer wall after long turns and remains near the outer wall during the straightaways. 

## Results

### Animation
![Vehicle animation](figures/animation.gif)

### Lap Statistics

| Metric | Value |
|---|---|
| Target speed | 15 m/s |
| Laps completed | 1 |
| Lap time | 206.0 s |
| Total distance traveled | 3 900 m |
| Average speed | 15.00 m/s |
| Off-track events | 0 — car stayed within boundaries |

The vehicle completed one full lap (3 057 m center-line length) in **206 seconds** at a constant 15 m/s. The additional ~843 m beyond the lap length reflects the continued simulation, and the car's continued trajectory past the finish line. The vehicle stayed inside the track boundaries for the entire simulation.

---

### Track Layout
![Track overview](figures/track_overview.png)

---

### Simulated Vehicle Path
The green line shows the simulated vehicle path over one full lap. Blue rectangles show vehicle pose at eight evenly-spaced time steps. The path stays within the track boundaries (blue) on all four sections.

![Vehicle path](figures/vehicle_path.png)

---

### Vehicle Heading vs. Time
Heading increases from 0° (east) to 180° (west) through the right curve (~t = 60–105 s), holds at 180° along the top straight (~t = 105–165 s), then increases to 360° through the left curve (~t = 165–210 s), and finally continues back along the main straightaway at 360°. 

![Heading over time](figures/heading.png)

---

### Actual Speed vs. Time
Tangential speed is held constant at the target 15 m/s by the driver model. Lateral dynamics do not affect the commanded forward speed (yet!).

![Speed over time](figures/speed.png)

## How to Run
1. Open MATLAB, `cd` into `Project2/`.
2. Run `run.m` — this calls `init`, `gentrack`, then `animate` (launches Simulink and animates).
3. Lap statistics (loops completed, lap times, track violations) are printed to the console at the end.
4. To regenerate figures: type `gen_figs_week2` from the `Project2/` directory.
5. To regenerate the animation GIF: type `gen_gif` from the `Project2/` directory.
