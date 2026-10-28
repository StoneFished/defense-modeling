classdef ScenarioTest < matlab.unittest.TestCase
    %SCENARIOTEST Truth trajectory: shapes, reproducibility and kinematics.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            srcDir = fullfile(fileparts(fileparts(mfilename("fullpath"))), "src");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));
        end
    end

    methods (Test)
        function outputShapes(testCase)
            [X, t, omega] = simulateScenario(0.5, 80, 1);
            testCase.verifySize(X, [6 80]);
            testCase.verifySize(omega, [1 80]);
            testCase.verifyEqual(t, (0:79)*0.5);
        end

        function seedIsReproducible(testCase)
            testCase.verifyEqual(simulateScenario(1, 50, 4), simulateScenario(1, 50, 4));
            testCase.verifyNotEqual(simulateScenario(1, 50, 4), simulateScenario(1, 50, 5));
        end

        function straightLegsAreConstantVelocity(testCase)
            T = 1;
            [X, ~, omega] = simulateScenario(T, 120, 1);
            k = find(omega(1:end-1) == 0);
            testCase.assertNotEmpty(k);
            testCase.verifyEqual(X(4:6, k+1), X(4:6, k), "AbsTol", 1e-9);
            testCase.verifyEqual(X(1:3, k+1), X(1:3, k) + T*X(4:6, k), "AbsTol", 1e-6);
        end

        function turnIsCircularArc(testCase)
            % During a coordinated turn the speed, the climb rate and the
            % centre of the turn circle stay fixed. An Euler step would make
            % the centre drift by metres per step.
            T = 1;
            [X, ~, omega] = simulateScenario(T, 120, 1);
            k = find(omega ~= 0);
            testCase.assertNotEmpty(k);
            w = omega(k(1));
            kk = [k, k(end) + 1];                 % every state on the arc
            v = X(4:6, kk);
            centre = X(1:2, kk) + [-v(2, :); v(1, :)] / w;
            testCase.verifyEqual(centre - centre(:, 1), zeros(2, numel(kk)), "AbsTol", 1e-6, ...
                "turn centre moved");
            testCase.verifyEqual(vecnorm(v(1:2, :)) - norm(v(1:2, 1)), zeros(1, numel(kk)), "AbsTol", 1e-9, ...
                "horizontal speed changed");
            testCase.verifyEqual(v(3, :) - v(3, 1), zeros(1, numel(kk)), "AbsTol", 1e-9, ...
                "climb rate changed");
            testCase.verifyEqual(X(3, k+1), X(3, k) + T*X(6, k), "AbsTol", 1e-6, ...
                "altitude does not follow the climb rate");
        end

        function turnChangesHeadingByTurnAngle(testCase)
            for a = [pi/2, -pi/3, 3*pi/4]
                X = simulateScenario(1, 90, 2, TurnAngle=a);
                headingChange = wrapAngle(atan2(X(5, end), X(4, end)) - atan2(X(5, 1), X(4, 1)));
                testCase.verifyEqual(headingChange, a, "AbsTol", 1e-9, sprintf("TurnAngle = %g", a));
            end
        end

        function zeroTurnAngleIsConstantVelocity(testCase)
            T = 2;
            [X, t] = simulateScenario(T, 40, 3, TurnAngle=0);
            testCase.verifyEqual(X(4:6, :) - X(4:6, 1), zeros(3, 40), "AbsTol", 1e-9);
            testCase.verifyEqual(X(1:3, :), X(1:3, 1) + X(4:6, 1)*t, "AbsTol", 1e-6);
        end
    end
end
