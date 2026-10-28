function [xEst, PEst, modeProb] = immFilter(z, T, R, x0, P0, radarPos, params)
%IMMFILTER Interacting multiple model (IMM) filter: constant velocity + coordinated turn.
%   [xEst, PEst, modeProb] = immFilter(z, T, R, x0, P0, radarPos) filters
%   the radar measurements z (3 x N, columns [range; azimuth; elevation],
%   one every T seconds) from a radar at radarPos with two mode-matched
%   EKFs that share each measurement:
%
%       mode 1, CV   nearly constant velocity (cvModel), state [p; v]
%       mode 2, CT   nearly coordinated turn (ctModel), state [p; v; omega]
%
%   It returns, after each measurement, the combined estimate of [p; v]
%   xEst (6 x N), its covariance PEst (6 x 6 x N) and the mode
%   probabilities modeProb (2 x N, rows CV and CT).
%
%   R, x0 and P0 are as for ekf: x0, P0 is the prior on [p; v] at the time
%   of z(:,1), the first step is an update only, and every later step is
%   predict, then update.
%
%   immFilter(..., params) sets the process noise of both modes, the
%   mode-transition matrix and the initial mode probabilities; see
%   immParams for the fields and the defaults.
%
%   Each step has four parts:
%     1. Mixing: each mode restarts from a blend of both modes' last
%        estimates, weighted by the probability that the target was in
%        that mode given that it is in this one now.
%     2. Mode-matched filtering: one EKF predict and update per mode.
%     3. Mode probabilities: reweighted by how well each mode predicted
%        the measurement (its likelihood).
%     4. Combination: the probability-weighted mixture of the two modes.
%
%   The CV mode has no turn-rate state. Where its estimate stands in for
%   the CT mode (the prior, and mixing into the CT mode), omega is set to
%   0 with a standard deviation of 10 deg/s. That is wide enough for the
%   CT mode to pick up a new turn quickly, and it is the convention the
%   toolbox's trackingIMM uses (see toolboxImm).
%
%   See also immParams, ekfUpdate, cvModel, ctModel, ekf, toolboxImm.

arguments
    z (3,:) double
    T (1,1) double {mustBePositive}
    R (3,3) double
    x0 (6,1) double
    P0 (6,6) double
    radarPos (3,1) double
    params (1,1) struct = immParams()
end

sigmaOmegaCv = deg2rad(10);

N = size(z, 2);
xEst = zeros(6, N);
PEst = zeros(6, 6, N);
modeProb = zeros(2, N);

% TODO: the starting mode-conditioned estimates. The CV mode starts from
% (x0, P0); the CT mode from x0 and P0 extended with omega = 0 and
% variance sigmaOmegaCv^2. The mode probabilities start at
% params.ModeProb0.

for k = 1:N
    if k > 1
        % TODO 1. Mixing. With TM = params.TransitionMatrix and the last
        %   mode probabilities mu: the predicted mode probabilities
        %   muPred(j) = sum_i TM(i,j)*mu(i), the mixing weights
        %   TM(i,j)*mu(i)/muPred(j), and each mode's starting estimate:
        %   the mean and covariance (including the spread of the means)
        %   of the weighted mixture of both modes' last estimates. Convert
        %   between the 6- and 7-state spaces as described above.
        % TODO 2a. Predict each mode: cvModel for CV, ctModel for CT.
    end

    % TODO 2b. Update both modes with z(:, k) using ekfUpdate, keeping the
    %   innovation nu and its covariance S.
    % TODO 3. Mode probabilities: muPred (at k = 1, the starting
    %   probabilities) times each mode's likelihood N(nu; 0, S),
    %   normalized. Work with log-likelihoods: radar likelihoods can be
    %   small enough to underflow.
    % TODO 4. Combination: the probability-weighted mixture of the two
    %   modes in [p; v] coordinates, stored in xEst(:, k), PEst(:, :, k)
    %   and modeProb(:, k).
    error("immFilter:notImplemented", "Implement the IMM steps in immFilter.m.")
end
end
