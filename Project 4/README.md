# Project 4 (Team) – Week 1


In Project 4, we combined lateral and longitudinal dynamics with an electric powertrain model to simulate a vehicle completing laps around a track. The objective for Week 2 was to maximize the number of laps completed in a 60-minute simulation while satisfying all project constraints.
---

## Run Instructions

```matlab
% Full pipeline:
p4_init        % initialize parameters and generate track
p4_runsim      % run simulation and generate plots
```

### Key Scripts

| Script        | Purpose                                                                        |
| ------------- | ------------------------------------------------------------------------------ |
| `p4_init.m`   | Initializes vehicle, battery, and motor parameters; sets simulation conditions |
| `p4_runsim.m` | Runs simulation, computes performance metrics, generates plots                 |
| `gentrack.m`  | Generates oval track geometry                                                  |
| `friction.m`  | Enforces braking constraints                                                   |


---

## Results

| Parameter             | Value             |
| --------------------- | ----------------- |
| Simulation time       | 3600 s (60 min)   |
| Target speed          | 20 mph (8.94 m/s) |
| Lookahead distance    | 7 m               |
| Laps completed        | 10.43             |
| Track length          | 3057 m/lap        |
| Initial SOC           | 80.00%            |
| Final SOC             | 78.85%            |
| SOC drop              | 1.15%             |
| Max cross-track error | 7.27 m            |
| Track limit           | ±7.50 m           |
| Week 2 valid          | YES               |


---

## Figures

### Fig 1 — Lap 1 Trajectory

![Lap 1 Path](figures/fig1_lap1_path.png)

### Fig 2 — Full Run Trajectory (All Laps, coloured by time)

![All Laps Path](figures/fig2_all_laps_path.png)

### Fig 3 — Vehicle Speed vs Time (Full Run)

![Speed vs Time](figures/fig3_speed.png)

### Fig 3b — Vehicle Speed vs Time (First 250 s)

![Speed vs Time Zoom](figures/fig3b_speed_zoom.png)

### Fig 4 — Battery SOC vs Time

![Battery SOC](figures/fig4_soc.png)

### Fig 5 — Lateral Position vs Lap Position (Cross-Track Error)

![Lateral Laps](figures/fig5_lateral_laps.png)

### Animation

![Track Animation](figures/animation.gif)

---

## Observations

- The vehicle successfully completed 10.43 laps, exceeding the minimum requirement of 5 laps.
- A target speed of 20 mph was selected to ensure stability and prevent the vehicle from leaving the track.
- The maximum cross-track error (7.27 m) remained within the allowable track boundary of ±7.5 m, though it operated close to the limit.
- Battery SOC remained well within constraints, dropping only 1.15% over the full simulation.
- The constant-speed strategy ensured smooth and stable tracking, but limits maximum achievable laps due to conservative cornering speed.

---
Conclusion

The final Week 2 model meets all project requirements and achieves a valid simulation. The results demonstrate that lateral tracking performance is the primary limiting factor, as increasing speed beyond 20 mph caused the vehicle to leave the track. Future improvements could include adaptive speed control (slower in curves, faster on straights) to further increase lap count while maintaining stability.
