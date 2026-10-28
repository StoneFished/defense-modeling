function params = immParams(opts)
%IMMPARAMS Tuning of the CV/CT IMM: process noise, mode switching, start.
%   params = immParams() returns the Day 2 tuning as a struct for
%   immFilter and toolboxImm. immParams(Name=Value) overrides fields:
%
%     SigmaACv          CV mode: random acceleration per axis [m/s^2]
%     SigmaACt          CT mode: random acceleration per axis [m/s^2]
%     SigmaOmegaDot     CT mode: random turn-rate change [rad/s^2]
%     TransitionMatrix  2x2; element (i,j) is the probability that the
%                       target switches from mode i to mode j between two
%                       measurements (1 = CV, 2 = CT); rows sum to 1
%     ModeProb0         mode probabilities [CV; CT] at the first
%                       measurement
%
%   See also immFilter, toolboxImm.

arguments
    opts.SigmaACv (1,1) double {mustBeNonnegative} = 1
    opts.SigmaACt (1,1) double {mustBeNonnegative} = 1
    opts.SigmaOmegaDot (1,1) double {mustBeNonnegative} = deg2rad(6)
    opts.TransitionMatrix (2,2) double {mustBeNonnegative} = [0.99 0.01; 0.05 0.95]
    opts.ModeProb0 (2,1) double {mustBeNonnegative} = [0.9; 0.1]
end

if any(abs(sum(opts.TransitionMatrix, 2) - 1) > 1e-12) || abs(sum(opts.ModeProb0) - 1) > 1e-12
    error("immParams:notStochastic", "TransitionMatrix rows and ModeProb0 must each sum to 1.")
end
params = opts;
end
