classdef EkfTest < matlab.unittest.TestCase
    %EKFTEST Motion model checks and the Day 1 EKF acceptance checks.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            srcDir = fullfile(fileparts(fileparts(mfilename("fullpath"))), "src");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));
        end
    end

    methods (Test)
        function cvTransitionMovesAtConstantVelocity(testCase)
            T = 0.5;
            F = cvModel(T, 3);
            x = [1; 2; 3; 40; -50; 6];
            testCase.assertSize(F, [6 6]);
            testCase.verifyEqual(F*x, [x(1:3) + T*x(4:6); x(4:6)], "AbsTol", 1e-12);
        end

        function cvProcessNoiseIsCovariance(testCase)
            [~, Q] = cvModel(0.5, 3);
            testCase.assertSize(Q, [6 6]);
            testCase.assertGreaterThan(norm(Q), 0);
            tol = 1e-12*norm(Q);
            testCase.verifyEqual(Q, Q', "AbsTol", tol, "Q must be symmetric");
            testCase.verifyGreaterThan(min(eig((Q + Q')/2)), -tol, "Q must be positive semidefinite");
            [~, Q2] = cvModel(0.5, 6);
            testCase.verifyEqual(Q2, 4*Q, "AbsTol", 4*tol, "Q must scale with sigmaA^2");
        end

        function outputsAreWellFormed(testCase)
            steps = 30;
            [X, z, F, Q, R] = maneuveringRun(steps, 1);
            P0 = blkdiag(100^2*eye(3), 20^2*eye(3));
            [xEst, PEst] = ekf(z, F, Q, R, X(:, 1), P0, [0; 0; 0]);
            testCase.verifySize(xEst, [6 steps]);
            testCase.assertSize(PEst, [6 6 steps]);
            asym = max(abs(PEst - pagetranspose(PEst)), [], "all");
            minEig = min(arrayfun(@(k) min(eig((PEst(:, :, k) + PEst(:, :, k)')/2)), 1:steps));
            testCase.verifyLessThan(asym, 1e-6, "covariances must stay symmetric");
            testCase.verifyGreaterThan(minEig, 0, "covariances must stay positive definite");
        end

        function noiseFreeCvTrackLocksOntoTruth(testCase)
            % Acceptance check. Perfect measurements of a constant-velocity
            % target and a prior that is off by hundreds of metres: after a
            % burn-in the estimate must match the truth to numerical
            % precision. The filter's R is small but not ~0: as R -> 0 the
            % velocity error of a nearly-constant-velocity filter flips sign
            % every step instead of dying out. With these values the error
            % shrinks several-fold per step.
            T = 1;
            steps = 40;
            burnIn = 15;
            radarPos = [0; 0; 0];
            X = simulateScenario(T, steps, 1, TurnAngle=0);
            z = generateMeasurements(X, radarPos, zeros(3), 1);
            [F, Q] = cvModel(T, 10);
            R = diag([1.5, 6e-5, 6e-5].^2);         % ~1.5 m in range and cross-range
            x0 = X(:, 1) + [400; -300; 150; 25; -20; 5];
            P0 = blkdiag(1e3^2*eye(3), 100^2*eye(3));

            xEst = ekf(z, F, Q, R, x0, P0, radarPos);

            after = burnIn+1:steps;
            testCase.verifyLessThan(vecnorm(xEst(1:3, after) - X(1:3, after)), 1e-6, ...
                "position error [m] after burn-in");
            testCase.verifyLessThan(vecnorm(xEst(4:6, after) - X(4:6, after)), 1e-6, ...
                "velocity error [m/s] after burn-in");
        end

        function azimuthInnovationIsWrapped(testCase)
            % Extra check (not on the Day 1 list): a target flying along the
            % radar's -x axis, where measured azimuths flip between +pi and
            % -pi. An unwrapped innovation looks like a 2*pi error.
            T = 1;
            steps = 60;
            radarPos = [0; 0; 0];
            x0 = [-25e3; 30; 4e3; 150; 0; 0];
            X = x0 + [T*x0(4:6); zeros(3, 1)] * (0:steps-1);
            R = diag([20, deg2rad(0.2), deg2rad(0.2)].^2);
            z = generateMeasurements(X, radarPos, R, 2);
            [F, Q] = cvModel(T, 1);
            P0 = blkdiag(100^2*eye(3), 20^2*eye(3));

            xEst = ekf(z, F, Q, R, X(:, 1), P0, radarPos);

            testCase.verifyLessThan(max(vecnorm(xEst(1:3, :) - X(1:3, :))), 500, ...
                "max position error [m] near azimuth +/-pi");
        end

        function rmseBeatsRawMeasurements(testCase)
            % Acceptance check on the maneuvering scenario, pooled over
            % several noise draws so one lucky draw cannot decide it. The
            % constant-velocity model lags in the turn, which caps the gain:
            % a correct filter lands near 0.68 here (~0.55 on the straight
            % legs alone).
            steps = 120;
            radarPos = [0; 0; 0];
            P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
            measSq = 0;
            ekfSq = 0;
            for seed = 1:10
                [X, z, F, Q, R] = maneuveringRun(steps, seed);
                x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];
                xEst = ekf(z, F, Q, R, x0, P0, radarPos);
                measSq = measSq + sum(vecnorm(meas2cart(z, radarPos) - X(1:3, :)).^2);
                ekfSq = ekfSq + sum(vecnorm(xEst(1:3, :) - X(1:3, :)).^2);
            end
            ratio = sqrt(ekfSq / measSq);
            testCase.verifyLessThan(ratio, 0.8, ...
                sprintf("EKF position RMSE / raw measurement RMSE = %.2f", ratio));
        end
    end
end

function [X, z, F, Q, R] = maneuveringRun(steps, noiseSeed)
% The Day 1 setup: default scenario (radar at the origin), the same sensor
% and filter settings as experiments/exp01_ekf_track.m.
T = 1;
R = diag([25, deg2rad(0.2), deg2rad(0.2)].^2);
X = simulateScenario(T, steps, 1);
z = generateMeasurements(X, [0; 0; 0], R, noiseSeed);
[F, Q] = cvModel(T, 10);
end
