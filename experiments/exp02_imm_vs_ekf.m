%EXP02_IMM_VS_EKF Day 2: an IMM recovers from a maneuver much faster than a CV-EKF.
%   Monte Carlo over 100 trials of one scenario: straight flight, a 4 g,
%   90-degree turn from t = 40 s to 50 s, then straight again. Every trial
%   draws fresh measurement noise and fresh 0.1 g random accelerations of
%   the target. Runs the single-model CV-EKF (Day 1), the CV/CT IMM and,
%   as a benchmark, the toolbox's trackingIMM set up the same way. Prints
%   the headline numbers and saves figures/exp02_imm_vs_ekf.png (error vs
%   time) and figures/exp02_diagnostics.png (mode probability and NEES).

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

%% Scenario, radar and trials
T = 1;                              % sample period [s]
steps = 100;
radarPos = [0; 0; 0];
R = diag([25, deg2rad(0.2), deg2rad(0.2)].^2);  % same radar as Day 1
turnWindow = [40 50];               % [s]
turnAngle = pi/2;                   % 90 degrees in 10 s at ~250 m/s: ~4 g
accelNoise = 1;                     % target's random accelerations [m/s^2]
nTrials = 100;

%% Filters
% The CV-EKF is the IMM's CV mode on its own. Both are tuned to the
% target's 1 m/s^2 random accelerations, i.e. for straight flight.
params = immParams();
[F, Q] = cvModel(T, params.SigmaACv);
P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
runBenchmark = ~isempty(which("trackingIMM"));

%% Monte Carlo
[errEkf, errImm, neesEkf, neesImm, muCt] = deal(zeros(nTrials, steps));
benchDiff = nan(1, nTrials);        % max |immFilter - trackingIMM| position gap per trial [m]
for trial = 1:nTrials
    [X, t] = simulateScenario(T, steps, 1, TurnAngle=turnAngle, TurnWindow=turnWindow, ...
        AccelNoise=accelNoise, NoiseSeed=trial);
    z = generateMeasurements(X, radarPos, R, trial);
    x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];

    [xEkf, PEkf] = ekf(z, F, Q, R, x0, P0, radarPos);
    [xImm, PImm, modeProb] = immFilter(z, T, R, x0, P0, radarPos, params);
    errEkf(trial, :) = vecnorm(xEkf(1:3, :) - X(1:3, :));
    errImm(trial, :) = vecnorm(xImm(1:3, :) - X(1:3, :));
    neesEkf(trial, :) = nees(X, xEkf, PEkf);
    neesImm(trial, :) = nees(X, xImm, PImm);
    muCt(trial, :) = modeProb(2, :);

    if runBenchmark
        xTbx = toolboxImm(z, T, R, x0, P0, radarPos, params);
        benchDiff(trial) = max(vecnorm(xImm(1:3, :) - xTbx(1:3, :)));
    end
end

%% Headline numbers
onset = turnWindow(1);
before = t >= 20 & t < onset;       % straight leg, after the start-up transient
after = t >= onset;
turn = t >= onset & t < turnWindow(2);
settled = t >= turnWindow(2) + 10;  % straight again, 10 s after the turn
meanEkf = mean(errEkf);
meanImm = mean(errImm);
peakEkf = max(meanEkf(after));
peakImm = max(meanImm(after));
recoverEkf = timeToRecover(t, meanEkf, onset, 1.25*mean(meanEkf(before)));
recoverImm = timeToRecover(t, meanImm, onset, 1.25*mean(meanImm(before)));
meanMuCt = mean(muCt);
switchDelay = t(find(after & meanMuCt > 0.5, 1)) - onset;
aneesEkf = mean(neesEkf);
aneesImm = mean(neesImm);
% 95% band for the average of nTrials NEES values of a 6-state filter:
% nTrials*ANEES is chi-square with 6*nTrials degrees of freedom.
neesBand = 2*gammaincinv([0.025 0.975], 6*nTrials/2)/nTrials;

fprintf("Monte Carlo: %d trials, 4 g turn from t = %g s to %g s\n", nTrials, turnWindow);
fprintf("Before the turn:  mean position error CV-EKF %.1f m, IMM %.1f m (IMM %+.0f%%)\n", ...
    mean(meanEkf(before)), mean(meanImm(before)), 100*(mean(meanImm(before))/mean(meanEkf(before)) - 1));
fprintf("Peak after onset: CV-EKF %.0f m, IMM %.0f m (%.0f%% lower)\n", ...
    peakEkf, peakImm, 100*(1 - peakImm/peakEkf));
fprintf("Time to recover:  CV-EKF %g s, IMM %g s (%.0f%% faster)\n", ...
    recoverEkf, recoverImm, 100*(1 - recoverImm/recoverEkf));
fprintf("Mode switch: mean CT probability %.2f before the turn, above 0.5 %g s after onset\n", ...
    mean(meanMuCt(before)), switchDelay);
fprintf("ANEES (expected 6, 95%% band %.2f-%.2f):\n", neesBand);
fprintf("  before the turn  CV-EKF %5.1f   IMM %5.1f\n", mean(aneesEkf(before)), mean(aneesImm(before)));
fprintf("  in the turn      CV-EKF %5.0f   IMM %5.1f\n", mean(aneesEkf(turn)), mean(aneesImm(turn)));
fprintf("  10 s after it    CV-EKF %5.1f   IMM %5.1f\n", mean(aneesEkf(settled)), mean(aneesImm(settled)));
if runBenchmark
    fprintf("Benchmark: immFilter vs toolbox trackingIMM, largest position gap over all trials: %.1e m\n", ...
        max(benchDiff));
else
    fprintf("Benchmark skipped: trackingIMM needs the Sensor Fusion and Tracking Toolbox.\n");
end

%% Figures
colors.ink = [11 11 11]/255;
colors.secondary = [82 81 78]/255;
colors.muted = [137 135 129]/255;
colors.ekf = [42 120 214]/255;      % categorical slot 1, the CV-EKF as on Day 1
colors.imm = [235 104 52]/255;      % categorical slot 2
figDir = fullfile(repoRoot, "figures");
if ~isfolder(figDir)
    mkdir(figDir);
end

% Headline: mean position error vs time.
fig = figure(Name="exp02: IMM vs CV-EKF", Color="w", Position=[100 100 1000 520]);
if exist("theme", "file")
    theme(fig, "light");            % keep the saved figure light under a dark desktop theme
end
ax = axes(fig);
hold(ax, "on")
markTurn(ax, turnWindow, colors);
hEkf = plot(ax, t, meanEkf, "-", Color=colors.ekf, LineWidth=2, DisplayName="CV-EKF (single model)");
hImm = plot(ax, t, meanImm, "-", Color=colors.imm, LineWidth=2, DisplayName="IMM (CV + CT)");
[~, kPeak] = max(meanEkf);
text(ax, t(kPeak) + 2, peakEkf, sprintf("CV-EKF peak %.0f m", peakEkf), Color=colors.ink, ...
    VerticalAlignment="bottom")
[~, kPeakImm] = max(meanImm .* after);
text(ax, t(kPeakImm) + 2, peakImm*1.15, sprintf("IMM peak %.0f m", peakImm), Color=colors.ink, ...
    VerticalAlignment="bottom")
hold(ax, "off")
set(ax, YScale="log", YTick=[20 50 100 200 500 1000 2000], YLim=[20 2000], XLim=[0 t(end)])
ytickformat(ax, "%g")
grid(ax, "on")
ax.YMinorGrid = "off";
xlabel(ax, "Time [s]")
ylabel(ax, "Mean position error [m] (log scale)")
title(ax, sprintf("IMM vs single-model EKF through a 4 g turn, %d Monte Carlo trials", nTrials), ...
    sprintf("After the turn starts, the IMM's peak error is %.0f%% lower and it is back to normal in %g s instead of %g s", ...
    100*(1 - peakImm/peakEkf), recoverImm, recoverEkf), FontWeight="normal");
ax.Title.FontWeight = "bold";
legend(ax, [hEkf, hImm], Location="northeast")
exportgraphics(fig, fullfile(figDir, "exp02_imm_vs_ekf.png"), Resolution=150);

% Sanity checks: the IMM's mode probability and both filters' NEES.
fig = figure(Name="exp02: diagnostics", Color="w", Position=[150 150 1000 700]);
if exist("theme", "file")
    theme(fig, "light");
end
layout = tiledlayout(fig, 2, 1, TileSpacing="compact", Padding="compact");

ax = nexttile(layout);
hold(ax, "on")
markTurn(ax, turnWindow, colors);
band = prctile10to90(muCt);
fill(ax, [t, fliplr(t)], [band(1, :), fliplr(band(2, :))], colors.imm, FaceAlpha=0.15, ...
    EdgeColor="none", DisplayName="10-90% of trials");
plot(ax, t, meanMuCt, "-", Color=colors.imm, LineWidth=2, DisplayName="Mean over trials");
hold(ax, "off")
set(ax, YLim=[0 1], XLim=[0 t(end)])
grid(ax, "on")
ylabel(ax, "IMM probability of the CT mode")
title(ax, "Mode probability: the IMM switches to the turn model within seconds of the onset")
legend(ax, Location="northeast")

ax = nexttile(layout);
hold(ax, "on")
markTurn(ax, turnWindow, colors);
yregion(ax, neesBand(1), neesBand(2), FaceColor=colors.muted, FaceAlpha=0.25, ...
    DisplayName="95% band for an honest covariance");
yline(ax, 6, "--", Color=colors.secondary, HandleVisibility="off");
plot(ax, t, aneesEkf, "-", Color=colors.ekf, LineWidth=2, DisplayName="CV-EKF");
plot(ax, t, aneesImm, "-", Color=colors.imm, LineWidth=2, DisplayName="IMM");
hold(ax, "off")
set(ax, YScale="log", YLim=[1 1e4], YTick=[1 6 10 100 1000 10000], XLim=[0 t(end)])
ytickformat(ax, "%g")
grid(ax, "on")
ax.YMinorGrid = "off";
xlabel(ax, "Time [s]")
ylabel(ax, "Average NEES (expected 6, log scale)")
title(ax, "Consistency: NEES near 6 except while a filter's motion model is wrong")
legend(ax, Location="northeast")
exportgraphics(fig, fullfile(figDir, "exp02_diagnostics.png"), Resolution=150);

function markTurn(ax, turnWindow, colors)
% Shade the turn and mark its onset.
xregion(ax, turnWindow(1), turnWindow(2), FaceColor=colors.muted, FaceAlpha=0.12, ...
    HandleVisibility="off");
xline(ax, turnWindow(1), "-", "Turn starts", Color=colors.secondary, LabelVerticalAlignment="top", ...
    LabelOrientation="horizontal", HandleVisibility="off");
end

function tr = timeToRecover(t, err, onset, threshold)
% Time from the onset until err drops below threshold for good.
kLast = find(t >= onset & err > threshold, 1, "last");
if isempty(kLast)
    tr = 0;
else
    tr = t(kLast) + (t(2) - t(1)) - onset;
end
end

function band = prctile10to90(x)
% 10th and 90th percentiles of each column (base MATLAB, no toolbox).
xs = sort(x, 1);
n = size(x, 1);
band = xs(max(1, round([0.1; 0.9]*n)), :);
end
