classdef ImmTest < matlab.unittest.TestCase
    %IMMTEST IMM structure checks and the Day 2 acceptance checks.

    properties
        Mc  % Monte Carlo results shared by the acceptance checks
    end

    methods (TestClassSetup)
        function setUpMonteCarlo(testCase)
            % Add src/ to the path, then run the Day 2 Monte Carlo once (50
            % trials of the CV-EKF and the IMM) for the acceptance checks.
            srcDir = fullfile(fileparts(fileparts(mfilename("fullpath"))), "src");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));

            nTrials = 50;
            radarPos = [0; 0; 0];
            params = immParams();
            [F, Q] = cvModel(1, params.SigmaACv);
            P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
            mc = struct();
            [mc.errEkf, mc.errImm, mc.neesEkf, mc.neesImm, mc.muCt] = deal(zeros(nTrials, 100));
            for trial = 1:nTrials
                [X, z, R, t] = day2Trial(trial);
                x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];
                [xEkf, PEkf] = ekf(z, F, Q, R, x0, P0, radarPos);
                [xImm, PImm, modeProb] = immFilter(z, 1, R, x0, P0, radarPos, params);
                mc.errEkf(trial, :) = vecnorm(xEkf(1:3, :) - X(1:3, :));
                mc.errImm(trial, :) = vecnorm(xImm(1:3, :) - X(1:3, :));
                mc.neesEkf(trial, :) = nees(X, xEkf, PEkf);
                mc.neesImm(trial, :) = nees(X, xImm, PImm);
                mc.muCt(trial, :) = modeProb(2, :);
            end
            mc.t = t;
            testCase.Mc = mc;
        end
    end

    methods (Test)
        function outputsAreWellFormed(testCase)
            [X, z, R] = day2Trial(1);
            P0 = blkdiag(100^2*eye(3), 20^2*eye(3));
            [xEst, PEst, modeProb] = immFilter(z, 1, R, X(:, 1), P0, [0; 0; 0]);
            n = size(z, 2);
            testCase.verifySize(xEst, [6 n]);
            testCase.assertSize(PEst, [6 6 n]);
            testCase.assertSize(modeProb, [2 n]);
            testCase.verifyEqual(sum(modeProb), ones(1, n), "AbsTol", 1e-12);
            testCase.verifyGreaterThanOrEqual(modeProb, 0);
            asym = max(abs(PEst - pagetranspose(PEst)), [], "all");
            minEig = min(arrayfun(@(k) min(eig((PEst(:, :, k) + PEst(:, :, k)')/2)), 1:n));
            testCase.verifyLessThan(asym, 1e-6, "covariances must stay symmetric");
            testCase.verifyGreaterThan(minEig, 0, "covariances must stay positive definite");
        end

        function reducesToCvEkfWithoutTurnMode(testCase)
            % With the CT mode switched off (never entered, zero prior
            % probability) the IMM is exactly the CV-EKF.
            radarPos = [0; 0; 0];
            [~, z, R] = day2Trial(2);
            x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];
            P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
            params = immParams(TransitionMatrix=eye(2), ModeProb0=[1; 0]);
            [F, Q] = cvModel(1, params.SigmaACv);
            [xEkf, PEkf] = ekf(z, F, Q, R, x0, P0, radarPos);
            [xImm, PImm, modeProb] = immFilter(z, 1, R, x0, P0, radarPos, params);
            testCase.verifyEqual(xImm, xEkf, "AbsTol", 1e-9);
            testCase.verifyEqual(PImm, PEkf, "AbsTol", 1e-9*max(abs(PEkf), [], "all"));
            testCase.verifyEqual(modeProb(2, :), zeros(1, size(z, 2)));
        end

        function noHarmBeforeTheManeuver(testCase)
            % Acceptance check: on the straight leg before the turn, adding
            % the CT mode must not cost accuracy. Pooled over the trials, the
            % IMM's position RMSE stays within 10% of the CV-EKF's (it is
            % about 3% higher: a little CT probability is always mixed in).
            before = testCase.Mc.t >= 20 & testCase.Mc.t < 40;
            rmse = @(e) sqrt(mean(e.^2, "all"));
            ratio = rmse(testCase.Mc.errImm(:, before)) / rmse(testCase.Mc.errEkf(:, before));
            testCase.verifyLessThan(ratio, 1.1, ...
                sprintf("IMM / CV-EKF position RMSE before the turn = %.3f", ratio));
        end

        function switchesToTurnModeAfterOnset(testCase)
            % Acceptance check: the mean CT probability stays low on the
            % straight leg and passes 0.5 within 6 s of the turn starting.
            t = testCase.Mc.t;
            meanMuCt = mean(testCase.Mc.muCt);
            testCase.verifyLessThan(mean(meanMuCt(t >= 20 & t < 40)), 0.1, ...
                "mean CT probability before the turn");
            switchDelay = t(find(t >= 40 & meanMuCt > 0.5, 1)) - 40;
            testCase.assertNotEmpty(switchDelay, "CT probability never passed 0.5");
            testCase.verifyLessThanOrEqual(switchDelay, 6, ...
                "seconds from turn onset until the CT probability passes 0.5");
        end

        function recoversFasterThanCvEkf(testCase)
            % The Day 2 headline: after the onset the IMM's mean error peaks
            % far lower (about 92% lower here), and 5-10 s after the turn ends
            % it is back near its straight-flight level while the CV-EKF's is
            % still several times its own.
            t = testCase.Mc.t;
            meanEkf = mean(testCase.Mc.errEkf);
            meanImm = mean(testCase.Mc.errImm);
            before = t >= 20 & t < 40;
            peakRatio = max(meanImm(t >= 40)) / max(meanEkf(t >= 40));
            testCase.verifyLessThan(peakRatio, 0.3, ...
                sprintf("IMM / CV-EKF peak error after onset = %.2f", peakRatio));
            late = t >= 55 & t < 60;
            testCase.verifyLessThan(max(meanImm(late)), 1.25*mean(meanImm(before)), ...
                "IMM error 5-10 s after the turn, vs 1.25x its straight-flight level");
            testCase.verifyGreaterThan(min(meanEkf(late)), 3*mean(meanEkf(before)), ...
                "CV-EKF error 5-10 s after the turn, vs 3x its straight-flight level");
        end

        function neesIsNearItsExpectedValue(testCase)
            % Acceptance check: with honest covariances the average NEES of
            % the 6-state estimate is 6. Away from the maneuver the CV-EKF
            % sits at ~5.6 and the IMM at ~5: slightly conservative, because
            % the CT mode's turn-rate uncertainty widens the mixture. In the
            % turn the CV-EKF's model is wrong and its NEES runs into the
            % thousands; the IMM's stays an order of magnitude lower.
            t = testCase.Mc.t;
            aneesEkf = mean(testCase.Mc.neesEkf);
            aneesImm = mean(testCase.Mc.neesImm);
            straight = (t >= 20 & t < 40) | t >= 60;
            testCase.verifyGreaterThan(mean(aneesImm(straight)), 4, "IMM ANEES on the straight legs");
            testCase.verifyLessThan(mean(aneesImm(straight)), 8, "IMM ANEES on the straight legs");
            testCase.verifyGreaterThan(mean(aneesEkf(t >= 20 & t < 40)), 4.5, "CV-EKF ANEES before the turn");
            testCase.verifyLessThan(mean(aneesEkf(t >= 20 & t < 40)), 7.5, "CV-EKF ANEES before the turn");
            turn = t >= 40 & t < 50;
            testCase.verifyLessThan(mean(aneesImm(turn)), mean(aneesEkf(turn))/10, ...
                "IMM ANEES in the turn vs the CV-EKF's");
        end

        function matchesToolboxTrackingImm(testCase)
            % Benchmark: the toolbox's trackingIMM, set up with the same
            % models (toolboxImm), must produce the same estimates.
            testCase.assumeNotEmpty(which("trackingIMM"), ...
                "needs the Sensor Fusion and Tracking Toolbox");
            radarPos = [0; 0; 0];
            P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
            for trial = 1:3
                [~, z, R] = day2Trial(trial);
                x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];
                [xImm, ~, muImm] = immFilter(z, 1, R, x0, P0, radarPos);
                [xTbx, ~, muTbx] = toolboxImm(z, 1, R, x0, P0, radarPos);
                testCase.verifyLessThan(max(vecnorm(xImm(1:3, :) - xTbx(1:3, :))), 1e-3, ...
                    sprintf("trial %d: largest position gap [m]", trial));
                testCase.verifyLessThan(max(abs(muImm - muTbx), [], "all"), 1e-6, ...
                    sprintf("trial %d: largest mode probability gap", trial));
            end
        end
    end
end

function [X, z, R, t] = day2Trial(seed)
% One Monte Carlo trial of the Day 2 setup (experiments/exp02_imm_vs_ekf.m):
% radar at the origin, a 4 g, 90-degree turn from t = 40 s to 50 s, 0.1 g
% random accelerations, fresh process and measurement noise per seed.
R = diag([25, deg2rad(0.2), deg2rad(0.2)].^2);
[X, t] = simulateScenario(1, 100, 1, TurnAngle=pi/2, TurnWindow=[40 50], AccelNoise=1, NoiseSeed=seed);
z = generateMeasurements(X, [0; 0; 0], R, seed);
end
