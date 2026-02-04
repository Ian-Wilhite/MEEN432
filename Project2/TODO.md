# Project 2 TODO

- [ ] Wire `init.m` outputs into the Simulink model (use the generated `carDataBus`) so vehicle parameters and initial conditions feed the subsystems instead of hard-coded values.
- [ ] In Simulink, complete the control loop: add/finish `Driver Model` to compute steering `delta_f`, `Lateral Dynamics` to propagate `v_y` and `psi` with tire slip/forces, and `Transformation/Rotation` to output X/Y in world coordinates.
- [ ] Add `To Workspace` blocks for X, Y, psi, velocity, and time so `animate.m` can read them; align signal names with the variables referenced in `animate.m` (`simout.X`, `simout.Y`, `simout.psi`, `simout.tout`).
- [ ] Implement path-following logic to hit the target speed while respecting the 15 m track width; tune gains/limits to maximize lap speed without leaving the track.
- [ ] Ensure the XY Graph in Simulink plots the path in real time and matches the generated waypoints from `gentrack.m`.
- [x] Export track waypoints (`path` struct) from MATLAB into Simulink (Data Store, From Workspace, or Bus Creator) so the model uses the same track definition as the animation script. (done via `assignin` in `gentrack.m`)
- [x] Integrate `raceStat` scoring after the Simulink run (can keep the helper inside `animate.m` or move to its own file) and fix the formatting of the loop print statement. (moved to `raceStat.m` with cleaned output)
- [x] Verify the animated patch rotates correctly with heading and that car dimensions/offsets make sense relative to track width; consider leaving motion trail for visualization. (trail added, patch retained)
- [x] Add error handling/guardrails in scripts (`run.m`, `animate.m`) for missing simulation data and empty Simulink runs. (try/catch + signal checks)
- [x] Document run instructions and assumptions in `README.md` (MATLAB/Simulink version, required toolboxes, how to change speeds/params) and include final performance notes when available. (README added; final performance to fill after tuning)
