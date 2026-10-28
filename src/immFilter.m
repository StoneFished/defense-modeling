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

[F, Qcv] = cvModel(T, params.SigmaACv);
TM = params.TransitionMatrix;
mu = params.ModeProb0;
xCv = x0;
PCv = P0;
xCt = [x0; 0];
PCt = blkdiag(P0, sigmaOmegaCv^2);

for k = 1:N
    if k > 1
        % 1. Mixing. muPred(j) is the probability of mode j before the
        % measurement; mixW(i,j) the probability of mode i at the last
        % step given mode j now. A mode that cannot be active (muPred = 0)
        % keeps its own estimate.
        muPred = TM'*mu;
        mixW = TM .* mu ./ muPred';
        I2 = eye(2);
        mixW(:, muPred == 0) = I2(:, muPred == 0);
        xCvAsCt = [xCv; 0];
        PCvAsCt = blkdiag(PCv, sigmaOmegaCv^2);
        [xCv0, PCv0] = mixGaussians([xCv, xCt(1:6)], cat(3, PCv, PCt(1:6, 1:6)), mixW(:, 1));
        [xCt0, PCt0] = mixGaussians([xCvAsCt, xCt], cat(3, PCvAsCt, PCt), mixW(:, 2));

        % 2a. Mode-matched prediction.
        xCv = F*xCv0;
        PCv = F*PCv0*F' + Qcv;
        [xCt, Fct, Qct] = ctModel(xCt0, T, params.SigmaACt, params.SigmaOmegaDot);
        PCt = Fct*PCt0*Fct' + Qct;
    else
        muPred = mu;
    end

    % 2b. Mode-matched update with the shared measurement.
    [xCv, PCv, nuCv, SCv] = ekfUpdate(xCv, PCv, z(:, k), R, radarPos);
    [xCt, PCt, nuCt, SCt] = ekfUpdate(xCt, PCt, z(:, k), R, radarPos);

    % 3. Mode probabilities, computed in logs so that tiny likelihoods
    % cannot underflow.
    logW = log(muPred) + [logGaussian(nuCv, SCv); logGaussian(nuCt, SCt)];
    mu = exp(logW - max(logW));
    mu = mu/sum(mu);

    % 4. Combination, in the CV mode's [p; v] coordinates.
    [xEst(:, k), PEst(:, :, k)] = mixGaussians([xCv, xCt(1:6)], cat(3, PCv, PCt(1:6, 1:6)), mu);
    modeProb(:, k) = mu;
end
end

function [x, P] = mixGaussians(xs, Ps, w)
% Mean and covariance of the Gaussian mixture sum_i w(i)*N(xs(:,i), Ps(:,:,i)).
x = xs*w;
P = zeros(size(Ps, 1));
for i = 1:numel(w)
    d = xs(:, i) - x;
    P = P + w(i)*(Ps(:, :, i) + d*d');
end
end

function l = logGaussian(nu, S)
% log N(nu; 0, S).
L = chol(S, "lower");
r = L \ nu;
l = -(r'*r)/2 - sum(log(diag(L))) - numel(nu)/2*log(2*pi);
end
