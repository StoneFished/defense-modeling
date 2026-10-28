function [F, Q] = cvModel(T, sigmaA)
%CVMODEL Discrete-time nearly-constant-velocity motion model.
%   [F, Q] = cvModel(T, sigmaA) returns the state transition matrix F (6x6)
%   over one sample period T [s] and the process noise covariance Q (6x6)
%   for the state x = [px; py; pz; vx; vy; vz], so that
%
%       x(k+1) = F*x(k) + w(k),   w(k) ~ N(0, Q).
%
%   The process noise stands in for unknown target accelerations with
%   standard deviation sigmaA [m/s^2], independent on each axis. sigmaA is
%   the knob that trades noise smoothing against lag during maneuvers.
%
%   See also ekf.

arguments
    T (1,1) double {mustBePositive}
    sigmaA (1,1) double {mustBeNonnegative}
end

% Discrete white-noise acceleration: the acceleration a is constant over
% each sample period and independent between periods and axes,
% a ~ N(0, sigmaA^2*I). Over one period it moves the state by G*a.
I = eye(3);
F = [I, T*I; zeros(3), I];
G = [T^2/2*I; T*I];
Q = sigmaA^2 * (G*G');
end
