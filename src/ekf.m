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
%   The measurement model is measModel, linearized by measJacobian.
%
%   See also cvModel, measModel, measJacobian, wrapAngle.

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
        % TODO: predict x and P forward one step.
    end

    % TODO: update x and P with the measurement z(:, k).
    %   - Linearize the measurement model at the predicted state.
    %   - Wrap the azimuth innovation with wrapAngle: a target near azimuth
    %     +/-pi otherwise produces a 2*pi "innovation".
    %   - The textbook covariance update can lose symmetry and positive
    %     definiteness numerically; consider the Joseph form.
    error("ekf:notImplemented", "Implement the predict/update steps in ekf.m.")

    xEst(:, k) = x;
    PEst(:, :, k) = P;
end
end
