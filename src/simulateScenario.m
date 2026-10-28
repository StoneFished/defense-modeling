function [X, t, omega] = simulateScenario(T, steps, seed, opts)
%SIMULATESCENARIO True trajectory of one target: straight, coordinated turn, straight.
%   [X, t, omega] = simulateScenario(T, steps, seed) returns the true state
%   history X (6 x steps) sampled every T seconds at times t (1 x steps,
%   starting at 0). The target flies straight for the first third of the
%   run, turns at a constant rate during the middle third, then flies
%   straight again. omega (1 x steps) is the turn rate [rad/s] applied from
%   t(k) to t(k+1).
%
%   simulateScenario(..., TurnAngle=a) sets the total heading change of the
%   turn [rad]; positive turns counter-clockwise seen from above. The
%   default is pi/2. TurnAngle=0 gives a constant-velocity trajectory.
%
%   simulateScenario(..., TurnWindow=[tStart tEnd]) turns from tStart up to
%   tEnd [s] instead of during the middle third of the run.
%
%   simulateScenario(..., AccelNoise=sigma) adds random accelerations with
%   standard deviation sigma [m/s^2] per axis, constant over each sample
%   period (the noise model of cvModel and ctModel), so the target flies
%   only nearly straight and turns only nearly coordinated. NoiseSeed=s
%   seeds them (default: seed); different noise seeds give different
%   realizations of the same scenario. The default sigma is 0.
%
%   State convention used throughout this repo (Cartesian, z up):
%       x = [px; py; pz; vx; vy; vz]   [m, m/s]
%
%   The seed jitters the initial position, speed and heading, so different
%   seeds give different geometries and the same seed gives the same truth.
%
%   See also generateMeasurements.

arguments
    T (1,1) double {mustBePositive}
    steps (1,1) double {mustBeInteger, mustBePositive}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.TurnAngle (1,1) double = pi/2
    opts.TurnWindow (1,2) double = [1 2] * (steps - 1) * T / 3
    opts.AccelNoise (1,1) double {mustBeNonnegative} = 0
    opts.NoiseSeed (1,1) double {mustBeInteger, mustBeNonnegative} = seed
end

% Nominal start: ~30 km from a radar at the origin, 5 km up, flying toward
% -x at ~250 m/s in a gentle descent.
rs = RandStream("mt19937ar", "Seed", seed);
u = 2*rand(rs, 5, 1) - 1;                       % uniform jitter in [-1, 1]
pos0 = [30e3; 10e3; 5e3] + [1e3; 1e3; 500] .* u(1:3);
speed = 250 + 20*u(4);
heading = pi + deg2rad(10)*u(5);
climbRate = -10;
x0 = [pos0; speed*cos(heading); speed*sin(heading); climbRate];

% Zero turn rate on the straight legs; a constant rate over the turn window
% (by default the middle third), scaled so the heading changes by exactly
% TurnAngle.
t = (0:steps-1) * T;
inTurn = t >= opts.TurnWindow(1) & t < opts.TurnWindow(2);
omega = zeros(1, steps);
omega(inTurn) = opts.TurnAngle / (T * max(nnz(inTurn), 1));

% Random accelerations, held constant over each step. They come from a
% different generator than the jitter above, so equal seeds still give
% independent draws.
ns = RandStream("mrg32k3a", "Seed", opts.NoiseSeed);
accel = opts.AccelNoise * randn(ns, 3, steps - 1);

X = zeros(6, steps);
X(:, 1) = x0;
for k = 1:steps-1
    X(:, k+1) = propagate(X(:, k), omega(k), T) + [T^2/2*accel(:, k); T*accel(:, k)];
end
end

function x = propagate(x, omega, T)
%PROPAGATE Advance the state x = [p; v] by T seconds at turn rate omega.
%   omega == 0: straight flight at constant velocity.
%   omega ~= 0: coordinated turn in the horizontal plane. The speed stays
%   constant, (vx, vy) rotates by omega*T (counter-clockwise for omega > 0)
%   and vz is unchanged. The update is exact (the rotating velocity is
%   integrated over the step), not an Euler step; tests/ScenarioTest.m
%   checks that the turn is a true circular arc.

p = x(1:3);
v = x(4:6);
if omega == 0
    p = p + v*T;
else
    % (vx, vy) rotates at rate omega. Integrating the rotating velocity over
    % [0, T] gives the horizontal displacement; the result divides by omega,
    % hence the separate straight-flight branch.
    s = sin(omega*T);
    c = cos(omega*T);
    p = p + [[s, c - 1; 1 - c, s]*v(1:2)/omega; v(3)*T];
    v(1:2) = [c, -s; s, c]*v(1:2);
end
x = [p; v];
end
