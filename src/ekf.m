function [xEst, PEst] = ekf(z, F, Q, R, x0, P0, radarPos)
%EKF Extended Kalman filter: linear motion model, nonlinear radar measurements.
%   [xEst, PEst] = ekf(z, F, Q, R, x0, P0, radarPos) filters the measurements
%   z (3 x N, columns [range; azimuth; elevation]) from a radar at radarPos
%   and returns the posterior state estimates xEst (6 x N) and covariances
%   PEst (6 x 6 x N), one per measurement.
%
%   F, Q    motion model x(k+1) = F*x(k) + w(k), w ~ N(0, Q) (see cvModel)
%   R       measurement noise covariance (3x3)
%   x0, P0  prior mean and covariance AT THE TIME OF z(:,1): the first step
%           is an update only; every later step is predict, then update.
%
%   Each update is ekfUpdate: measModel, linearized by measJacobian at the
%   predicted state.
%
%   See also cvModel, ekfUpdate, measModel, measJacobian, immFilter.

arguments
    z (3,:) double
    F (6,6) double
    Q (6,6) double
    R (3,3) double
    x0 (6,1) double
    P0 (6,6) double
    radarPos (3,1) double
end

N = size(z, 2);
xEst = zeros(6, N);
PEst = zeros(6, 6, N);
x = x0;
P = P0;

for k = 1:N
    if k > 1
        x = F*x;
        P = F*P*F' + Q;
    end
    [x, P] = ekfUpdate(x, P, z(:, k), R, radarPos);
    xEst(:, k) = x;
    PEst(:, :, k) = P;
end
end
