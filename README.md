# Radar Target Tracking for Missile Defense

A ground radar sees an incoming air target only as noisy range, azimuth and elevation measurements, and the target can maneuver at any moment. A defense system still has to know where the target is and where it is heading, accurately and in real time. This project builds that tracking pipeline step by step in MATLAB: it simulates target trajectories and radar measurements, estimates the target's 3D position and velocity with nonlinear Kalman filtering, and checks each piece against known answers. Day 1 is the baseline: an extended Kalman filter (EKF) with a nearly-constant-velocity motion model, tracking a target that flies straight, makes a coordinated turn, then flies straight again.

![Day 1: EKF tracking a maneuvering target from radar measurements](figures/exp01_ekf_track.png)

*Top: the scenario, and a zoom on the turn. The EKF estimate (blue) follows the truth (black) and smooths out the scatter of the raw measurements (gray). Bottom: position error over time. The constant-velocity model lags behind the target during the turn (shaded).*

## How to run

Developed on MATLAB R2026b; it should run on R2023a or newer. Day 1 uses base MATLAB only, no toolboxes.

From the repository root:

```matlab
runtests("tests")                      % unit and acceptance tests
run("experiments/exp01_ekf_track.m")   % simulate, track, save figures/exp01_ekf_track.png
```

## What the tests check

- The analytic measurement Jacobian matches central finite differences to a relative 1e-6.
- With noise-free measurements of a constant-velocity target, the EKF locks onto the truth to numerical precision after a short burn-in.
- On the maneuvering scenario, the EKF position RMSE is well below the error of the raw measurements.
- The coordinated turn is an exact circular arc, the measurement angles follow the documented conventions, and azimuth innovations are wrapped across ±π.

## Layout

| Path | Contents |
|---|---|
| `src/` | scenario simulation, radar measurement model and Jacobian, motion model, EKF |
| `experiments/` | one script per experiment; each saves its figure to `figures/` |
| `tests/` | `matlab.unittest` test classes |

Conventions: the state is `x = [px; py; pz; vx; vy; vz]` (m, m/s) in a Cartesian frame with z up. A measurement is `[range; azimuth; elevation]`, with azimuth measured in the x-y plane from +x toward +y and elevation measured up from that plane.
