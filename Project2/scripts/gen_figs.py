"""Generate track and vehicle-poses figures matching gentrack.m logic.
Run with: python3 scripts/gen_figs.py (from the Project2 directory).
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

radius = 200.0
l_st = 900.0
width = 15.0

delta_s = 10.0
l_curve = np.pi * radius
total_length = 2 * l_st + 2 * l_curve
npts = round(total_length / delta_s)
delta_s = total_length / (npts - 1)

delta_theta = delta_s / radius

xpath = np.zeros(npts)
ypath = np.zeros(npts)
tpath = np.zeros(npts)
xinpath = np.zeros(npts)
yinpath = np.zeros(npts)
xoutpath = np.zeros(npts)
youtpath = np.zeros(npts)

i = 0
while i < npts - 1:
    yinpath[0] = width / 2
    youtpath[0] = -width / 2
    if xpath[i] < l_st:
        if xpath[i] >= 0:
            if ypath[i] < radius:
                xpath[i + 1] = xpath[i] + delta_s
                ypath[i + 1] = ypath[i]
                xinpath[i + 1] = xinpath[i] + delta_s
                yinpath[i + 1] = yinpath[i]
                xoutpath[i + 1] = xoutpath[i] + delta_s
                youtpath[i + 1] = youtpath[i]
                tpath[i + 1] = 0
            else:
                xpath[i + 1] = xpath[i] - delta_s
                ypath[i + 1] = ypath[i]
                xinpath[i + 1] = xinpath[i] - delta_s
                yinpath[i + 1] = yinpath[i]
                xoutpath[i + 1] = xoutpath[i] - delta_s
                youtpath[i + 1] = youtpath[i]
                tpath[i + 1] = np.pi
        else:
            cx, cy = 0, radius
            rot = np.array([[np.cos(delta_theta), -np.sin(delta_theta)],
                            [np.sin(delta_theta), np.cos(delta_theta)]])
            rx, ry = rot @ np.array([xpath[i] - cx, ypath[i] - cy])
            rxin, ryin = rot @ np.array([xinpath[i] - cx, yinpath[i] - cy])
            rxout, ryout = rot @ np.array([xoutpath[i] - cx, youtpath[i] - cy])
            xpath[i + 1], ypath[i + 1] = rx + cx, ry + cy
            xinpath[i + 1], yinpath[i + 1] = rxin + cx, ryin + cy
            xoutpath[i + 1], youtpath[i + 1] = rxout + cx, ryout + cy
            tpath[i + 1] = tpath[i] + delta_theta
    else:
        cx, cy = l_st, radius
        rot = np.array([[np.cos(delta_theta), -np.sin(delta_theta)],
                        [np.sin(delta_theta), np.cos(delta_theta)]])
        rx, ry = rot @ np.array([xpath[i] - cx, ypath[i] - cy])
        rxin, ryin = rot @ np.array([xinpath[i] - cx, yinpath[i] - cy])
        rxout, ryout = rot @ np.array([xoutpath[i] - cx, youtpath[i] - cy])
        xpath[i + 1], ypath[i + 1] = rx + cx, ry + cy
        xinpath[i + 1], yinpath[i + 1] = rxin + cx, ryin + cy
        xoutpath[i + 1], youtpath[i + 1] = rxout + cx, ryout + cy
        tpath[i + 1] = tpath[i] + delta_theta
    i += 1

fig_dir = Path('figures')
fig_dir.mkdir(parents=True, exist_ok=True)

plt.figure(figsize=(8, 4))
plt.plot(xpath, ypath, '--r', label='Center line')
plt.plot(xinpath, yinpath, 'b', label='Inner border')
plt.plot(xoutpath, youtpath, 'b', label='Outer border')
plt.axis('equal')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Project 2 Track: 900 m straights, 200 m radius, 15 m width')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(fig_dir / 'track_overview.png', dpi=200)
plt.close()

L = 15
W = 5
sample_idx = np.linspace(0, len(xpath) - 1, 8, dtype=int)
plt.figure(figsize=(8, 4))
plt.plot(xpath, ypath, '--k', alpha=0.4)
plt.axis('equal')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Sample vehicle poses along center line')
plt.grid(True)
for idx in sample_idx:
    x, y = xpath[idx], ypath[idx]
    nxt = (idx + 1) % len(xpath)
    psi = np.arctan2(ypath[nxt] - y, xpath[nxt] - x)
    car = np.array([[-L / 2, -W / 2],
                    [-L / 2, W / 2],
                    [L / 2, W / 2],
                    [L / 2, -W / 2],
                    [-L / 2, -W / 2]])
    R = np.array([[np.cos(psi), -np.sin(psi)],
                  [np.sin(psi), np.cos(psi)]])
    rcar = car @ R.T + np.array([x, y])
    plt.plot(rcar[:, 0], rcar[:, 1], 'r-')
plt.tight_layout()
plt.savefig(fig_dir / 'vehicle_poses.png', dpi=200)
plt.close()

print('Figures written to', fig_dir.resolve())
