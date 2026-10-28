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

v = x(4:6);
omega = x(7);
s = sin(omega*T);
c = cos(omega*T);

% sw = sin(omega*T)/omega and cw = (1 - cos(omega*T))/omega integrate the
% rotating velocity; dsw and dcw are their derivatives with respect to
% omega. All four are smooth through omega = 0, where the direct formulas
% divide 0 by 0, so use their Taylor series near it.
wT = omega*T;
if abs(wT) < 1e-3
    sw = T*(1 - wT^2/6);
    cw = T*(wT/2 - wT^3/24);
    dsw = T^2*(-wT/3 + wT^3/30);
    dcw = T^2*(1/2 - wT^2/8);
else
    sw = s/omega;
    cw = (1 - c)/omega;
    dsw = (T*c - sw)/omega;
    dcw = (T*s - cw)/omega;
end

xPred = [x(1) + sw*v(1) - cw*v(2)
         x(2) + cw*v(1) + sw*v(2)
         x(3) + T*v(3)
         c*v(1) - s*v(2)
         s*v(1) + c*v(2)
         v(3)
         omega];

F = eye(7);
F(1, [4 5 7]) = [sw, -cw, dsw*v(1) - dcw*v(2)];
F(2, [4 5 7]) = [cw, sw, dcw*v(1) + dsw*v(2)];
F(3, 6) = T;
F(4, [4 5 7]) = [c, -s, -T*(s*v(1) + c*v(2))];
F(5, [4 5 7]) = [s, c, T*(c*v(1) - s*v(2))];

[~, Qcv] = cvModel(T, sigmaA);
Q = blkdiag(Qcv, (T*sigmaOmegaDot)^2);
end
