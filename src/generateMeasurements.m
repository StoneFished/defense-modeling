function z = generateMeasurements(X, radarPos, R, seed)
%GENERATEMEASUREMENTS Noisy radar measurements of a true trajectory.
%   z = generateMeasurements(X, radarPos, R, seed) applies measModel to each
%   column of the true states X (6 x N) and adds zero-mean Gaussian noise
%   with covariance R (3x3, units [m^2, rad^2, rad^2]). Returns z (3 x N),
%   rows [range; azimuth; elevation], with azimuth wrapped to [-pi, pi).
%   R = zeros(3) gives noise-free measurements. X itself is untouched, so
%   truth and measurements can be plotted side by side.
%
%   See also measModel, meas2cart.

arguments
    X (6,:) double
    radarPos (3,1) double
    R (3,3) double
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
end

N = size(X, 2);
z = zeros(3, N);
for k = 1:N
    z(:, k) = measModel(X(:, k), radarPos);
end

% Draw N(0, R) through a symmetric square root, which also works for a
% singular R such as zeros(3).
[V, D] = eig((R + R')/2);
L = V * diag(sqrt(max(diag(D), 0)));
rs = RandStream("mt19937ar", "Seed", seed);
z = z + L*randn(rs, 3, N);
z(2, :) = wrapAngle(z(2, :));
end
