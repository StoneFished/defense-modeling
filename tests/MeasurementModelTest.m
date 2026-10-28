classdef MeasurementModelTest < matlab.unittest.TestCase
    %MEASUREMENTMODELTEST Conventions of measModel and the Jacobian check.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            srcDir = fullfile(fileparts(fileparts(mfilename("fullpath"))), "src");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));
        end
    end

    methods (Test)
        function knownGeometry(testCase)
            radarPos = [100; -200; 50];
            vel = [30; -40; 5];   % h must not depend on velocity
            % target offset from radar   expected [range; azimuth; elevation]
            cases = {
                [1000; 0; 0]                [1000; 0; 0]
                [0; 1000; 0]                [1000; pi/2; 0]
                [0; -1000; 0]               [1000; -pi/2; 0]
                [1000; 1000; sqrt(2)*1000]  [2000; pi/4; pi/4]
                [0; 1000; -1000]            [sqrt(2)*1000; pi/2; -pi/4]
                };
            for i = 1:size(cases, 1)
                z = measModel([radarPos + cases{i, 1}; vel], radarPos);
                testCase.verifyEqual(z, cases{i, 2}, "AbsTol", 1e-9, ...
                    sprintf("target at radar + [%g %g %g]", cases{i, 1}));
            end

            % Straight down the -x axis the azimuth sits on the +/-pi cut.
            z = measModel([radarPos + [-1000; 0; 0]; vel], radarPos);
            testCase.verifyEqual(abs(z(2)), pi, "AbsTol", 1e-12);
        end

        function measModelRoundTripsThroughMeas2cart(testCase)
            radarPos = [100; -200; 50];
            rs = RandStream("mt19937ar", "Seed", 7);
            for i = 1:20
                p = radarPos + [60e3*(rand(rs, 2, 1) - 0.5); 15e3*rand(rs)];
                z = measModel([p; 0; 0; 0], radarPos);
                testCase.verifyEqual(meas2cart(z, radarPos), p, "AbsTol", 1e-6);
            end
        end

        function jacobianMatchesFiniteDifferences(testCase)
            % Acceptance check. Each row is compared relative to its own size:
            % the angle rows scale like 1/range, so an absolute 1e-6 would
            % hide real errors there.
            radarPos = [100; -200; 50];
            rs = RandStream("mt19937ar", "Seed", 3);
            for i = 1:25
                zPos = [2e3 + 48e3*rand(rs); pi*(2*rand(rs) - 1); -0.2 + 1.4*rand(rs)];
                x = [meas2cart(zPos, radarPos); 300*(rand(rs, 3, 1) - 0.5)];
                H = measJacobian(x, radarPos);
                testCase.assertSize(H, [3 6]);
                Hfd = finiteDifferenceJacobian(x, radarPos);
                relErr = max(abs(H - Hfd), [], 2) ./ vecnorm(Hfd, 2, 2);
                testCase.verifyLessThan(relErr, 1e-6, ...
                    sprintf("x = [%s]", join(string(round(x')), " ")));
            end
        end
    end
end

function J = finiteDifferenceJacobian(x, radarPos)
% Central differences of measModel. The azimuth difference is wrapped so a
% step across the +/-pi cut does not show up as a 2*pi jump.
h = 1e-2;   % [m] for position, [m/s] for velocity
J = zeros(3, 6);
for j = 1:6
    dx = zeros(6, 1);
    dx(j) = h;
    dz = measModel(x + dx, radarPos) - measModel(x - dx, radarPos);
    dz(2) = wrapAngle(dz(2));
    J(:, j) = dz / (2*h);
end
end
