function [x, P, nu, S] = ekfUpdate(x, P, z, R, radarPos)
%EKFUPDATE One EKF measurement update with a radar measurement.
%   [x, P] = ekfUpdate(x, P, z, R, radarPos) updates the state estimate x
%   and covariance P with one measurement z = [range; azimuth; elevation]
%   from a radar at radarPos with measurement noise covariance R. The
%   measurement model is measModel, linearized by measJacobian at x.
%
%   x must start with [px; py; pz; vx; vy; vz]. Further states, such as
%   ctModel's turn rate, do not enter the measurement; they are corrected
%   through their correlation with the kinematic states.
%
%   [x, P, nu, S] = ekfUpdate(...) also returns the innovation nu (azimuth
%   wrapped to [-pi, pi)) and its covariance S, from which immFilter
%   computes each mode's measurement likelihood.
%
%   See also ekf, immFilter, measModel, measJacobian.

arguments
    x (:,1) double
    P (:,:) double
    z (3,1) double
    R (3,3) double
    radarPos (3,1) double
end

n = numel(x);
H = [measJacobian(x(1:6), radarPos), zeros(3, n - 6)];

% Near azimuth +/-pi the raw azimuth difference can be off by 2*pi, so wrap it.
nu = z - measModel(x(1:6), radarPos);
nu(2) = wrapAngle(nu(2));
S = H*P*H' + R;
K = P*H'/S;
x = x + K*nu;

% Joseph form keeps P positive definite under rounding; averaging with
% its transpose removes round-off asymmetry.
A = eye(n) - K*H;
P = A*P*A' + K*R*K';
P = (P + P')/2;
end
