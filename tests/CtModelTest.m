classdef CtModelTest < matlab.unittest.TestCase
    %CTMODELTEST Coordinated-turn motion model: kinematics, Jacobian, noise.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            srcDir = fullfile(fileparts(fileparts(mfilename("fullpath"))), "src");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));
        end
    end

    methods (Test)
        function matchesTruthPropagation(testCase)
            % With no noise, ctModel must reproduce simulateScenario's path
            % step by step, through a 4 g turn and along the straight legs,
            % with the scenario's turn rate as the seventh state.
            T = 1;
            [X, ~, omega] = simulateScenario(T, 100, 1, TurnWindow=[40 50]);
            for k = 1:99
                xPred = ctModel([X(:, k); omega(k)], T, 0, 0);
                testCase.verifyEqual(xPred, [X(:, k+1); omega(k)], "AbsTol", 1e-6, ...
                    sprintf("step %d", k));
            end
        end

        function zeroTurnRateIsConstantVelocity(testCase)
            T = 0.5;
            x = [1e3; -2e3; 3e3; 40; -50; 6];
            [xPred, F, Q] = ctModel([x; 0], T, 3, 0.1);
            [Fcv, Qcv] = cvModel(T, 3);
            testCase.verifyEqual(xPred, [Fcv*x; 0], "AbsTol", 1e-9);
            testCase.verifyEqual(F(1:6, 1:6), Fcv, "AbsTol", 1e-12);
            testCase.verifyEqual(Q(1:6, 1:6), Qcv, "AbsTol", 1e-12);
        end

        function processNoiseIsCvBlockPlusTurnRate(testCase)
            T = 0.5;
            [~, ~, Q] = ctModel([1; 2; 3; 40; -50; 6; 0.1], T, 3, 0.02);
            [~, Qcv] = cvModel(T, 3);
            testCase.verifyEqual(Q, blkdiag(Qcv, (0.02*T)^2), "AbsTol", 1e-12);
        end

        function jacobianMatchesFiniteDifferences(testCase)
            % Including turn rates at and near zero, where ctModel switches to
            % Taylor series. Each row is compared relative to its own size.
            T = 1;
            rs = RandStream("mt19937ar", "Seed", 5);
            for omega = [0, 1e-6, 5e-4, 2e-3, 0.05, -0.2, 0.6]
                x = [40e3*(rand(rs, 3, 1) - 0.5); 300*(rand(rs, 3, 1) - 0.5); omega];
                [~, F] = ctModel(x, T, 1, 0.01);
                Ffd = finiteDifferenceJacobian(x, T);
                relErr = max(abs(F - Ffd), [], 2) ./ vecnorm(Ffd, 2, 2);
                testCase.verifyLessThan(relErr, 1e-6, sprintf("omega = %g", omega));
            end
        end
    end
end

function J = finiteDifferenceJacobian(x, T)
% Central differences of ctModel's propagation.
h = [1e-2*ones(6, 1); 1e-5];   % [m], [m/s], [rad/s]
J = zeros(7);
for j = 1:7
    dx = zeros(7, 1);
    dx(j) = h(j);
    J(:, j) = (ctModel(x + dx, T, 0, 0) - ctModel(x - dx, T, 0, 0)) / (2*h(j));
end
end
