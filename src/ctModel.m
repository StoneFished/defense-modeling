function [xPred, F, Q] = ctModel(x, T, sigmaA, sigmaOmegaDot)
%CTMODEL Discrete-time nearly-coordinated-turn motion model, turn rate as a state.
%   [xPred, F, Q] = ctModel(x, T, sigmaA, sigmaOmegaDot) propagates the
%   state x = [px; py; pz; vx; vy; vz; omega] by one sample period T [s].
%   It also returns the Jacobian F = d(xPred)/dx (7x7) at x, which the EKF
%   uses to propagate the covariance, and the process noise covariance Q
%   (7x7), so that
%
%       x(k+1) = f(x(k)) + w(k),   w(k) ~ N(0, Q).
%
%   The target turns in the horizontal plane at the rate omega [rad/s,
%   positive counter-clockwise seen from above] at constant speed, and
%   climbs or descends at constant vertical velocity. f is the same exact
%   closed form that simulateScenario uses for the truth, with omega now
%   part of the state.
%
%   Q is cvModel's random-acceleration noise (standard deviation sigmaA
%   [m/s^2] per axis) for [p; v], plus a random turn-rate change with
%   standard deviation sigmaOmegaDot [rad/s^2] for omega.
%
%   See also cvModel, immFilter, simulateScenario.

arguments
    x (7,1) double
    T (1,1) double {mustBePositive}
    sigmaA (1,1) double {mustBeNonnegative}
    sigmaOmegaDot (1,1) double {mustBeNonnegative}
end

% TODO: propagate x with the exact coordinated-turn closed form: the same
% kinematics as propagate() in simulateScenario.m, with omega taken from
% x(7) and carried over unchanged.
%   - Near omega = 0 the closed form divides 0 by 0. Use Taylor series
%     there; tests/CtModelTest.m checks omega = 0 and tiny turn rates.
%   - F is the Jacobian of the propagation. Its last column (how the
%     predicted state depends on omega) is the new part; derive it by hand.
%   - Q: cvModel's random-acceleration block for [p; v], plus the variance
%     of the turn-rate change over one period T.
error("ctModel:notImplemented", "Implement ctModel.m.")
end
