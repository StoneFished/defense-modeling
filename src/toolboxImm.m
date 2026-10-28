function [xEst, PEst, modeProb] = toolboxImm(z, T, R, x0, P0, radarPos, params)
%TOOLBOXIMM Benchmark: the immFilter setup run through the toolbox's trackingIMM.
%   [xEst, PEst, modeProb] = toolboxImm(z, T, R, x0, P0, radarPos, params)
%   takes the same inputs and returns the same outputs as immFilter, but
%   filters with the Sensor Fusion and Tracking Toolbox's trackingIMM: two
%   trackingEKF filters, one with the constvel and one with the constturn
%   motion model, both measuring range, azimuth and elevation through
%   cvmeas/ctmeas in a spherical frame. The motion models, process noise,
%   measurement noise, transition matrix and prior all match immFilter.
%
%   The toolbox orders the state [x; vx; y; vy; z; vz] (constturn puts the
%   turn rate, in deg/s, fifth) and its measurements [azimuth; elevation;
%   range] with angles in degrees; this function converts in and out of
%   the repo's conventions.
%
%   Requires the Sensor Fusion and Tracking Toolbox.
%
%   See also immFilter, immParams.

arguments
    z (3,:) double
    T (1,1) double {mustBePositive}
    R (3,3) double
    x0 (6,1) double
    P0 (6,6) double
    radarPos (3,1) double
    params (1,1) struct = immParams()
end

toTbx = [1 4 2 5 3 6];                  % repo [p; v] -> toolbox [x; vx; y; vy; z; vz]
fromTbx = [1 3 5 2 4 6];                % and back
D = diag([180/pi, 180/pi, 1]);          % [azimuth; elevation; range], degrees
zTbx = D*z([2 3 1], :);
RTbx = D*R([2 3 1], [2 3 1])*D;
measParams = struct("Frame", "spherical", "OriginPosition", radarPos, ...
    "Orientation", eye(3), "HasVelocity", false, "HasElevation", true);

% The process noise enters as accelerations [ax; ay; az] (constvel) and
% [ax; ay; omegaDot; az] (constturn, omegaDot in deg/s^2) with the same
% T^2/2 and T gains as cvModel and ctModel.
sa = params.SigmaACt;
cv = trackingEKF(@constvel, @cvmeas, zeros(6, 1), ...
    StateTransitionJacobianFcn=@constveljac, MeasurementJacobianFcn=@cvmeasjac, ...
    HasAdditiveProcessNoise=false, ProcessNoise=params.SigmaACv^2*eye(3), ...
    MeasurementNoise=RTbx, HasMeasurementWrapping=true);
ct = trackingEKF(@constturn, @ctmeas, zeros(7, 1), ...
    StateTransitionJacobianFcn=@constturnjac, MeasurementJacobianFcn=@ctmeasjac, ...
    HasAdditiveProcessNoise=false, ...
    ProcessNoise=diag([sa^2, sa^2, rad2deg(params.SigmaOmegaDot)^2, sa^2]), ...
    MeasurementNoise=RTbx, HasMeasurementWrapping=true);

% trackingIMM treats TransitionProbabilities as per second and raises them
% to the power dt; the T-th root makes one step match immFilter.
imm = trackingIMM({cv, ct}, TransitionProbabilities=params.TransitionMatrix^(1/T), ...
    ModelProbabilities=params.ModeProb0);

% The prior is given in the first (CV) filter's terms; switchimm converts
% it for the CT filter with omega = 0 and variance 100 (deg/s)^2.
initialize(imm, x0(toTbx), P0(toTbx, toTbx));

N = size(z, 2);
xEst = zeros(6, N);
PEst = zeros(6, 6, N);
modeProb = zeros(2, N);
for k = 1:N
    if k > 1
        predict(imm, T);
    end
    correct(imm, zTbx(:, k), measParams);
    xEst(:, k) = imm.State(fromTbx);
    PEst(:, :, k) = imm.StateCovariance(fromTbx, fromTbx);
    modeProb(:, k) = imm.ModelProbabilities;
end
end
