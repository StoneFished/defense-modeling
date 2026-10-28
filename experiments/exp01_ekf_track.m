%EXP01_EKF_TRACK Day 1: EKF tracking of one maneuvering target from one radar.
%   Simulates a target that flies straight, makes a coordinated turn and
%   flies straight again, generates noisy range/azimuth/elevation radar
%   measurements, runs a nearly-constant-velocity EKF, prints the position
%   RMSE and saves figures/exp01_ekf_track.png.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

%% Scenario and radar
T = 1;                              % sample period [s]
steps = 120;                        % number of radar measurements
radarPos = [0; 0; 0];
sigmaRange = 25;                    % [m]
sigmaAngle = deg2rad(0.2);          % azimuth and elevation [rad]
R = diag([sigmaRange, sigmaAngle, sigmaAngle].^2);

[X, t, omega] = simulateScenario(T, steps, 1);
z = generateMeasurements(X, radarPos, R, 1);

%% Filter
sigmaA = 10;                        % process noise: random acceleration [m/s^2]
[F, Q] = cvModel(T, sigmaA);

% Start at the first measurement with unknown velocity. P0 is wide, so the
% first few measurements quickly outweigh the prior.
x0 = [meas2cart(z(:, 1), radarPos); 0; 0; 0];
P0 = blkdiag(500^2*eye(3), 300^2*eye(3));
xEst = ekf(z, F, Q, R, x0, P0, radarPos);

%% Position errors
pMeas = meas2cart(z, radarPos);
measErr = vecnorm(pMeas - X(1:3, :));
ekfErr = vecnorm(xEst(1:3, :) - X(1:3, :));
rmse = @(e) sqrt(mean(e.^2));
fprintf("Position RMSE, raw measurements: %6.1f m\n", rmse(measErr));
fprintf("Position RMSE, EKF estimate:     %6.1f m\n", rmse(ekfErr));

%% Figure
colors.ink = [11 11 11]/255;        % truth and radar
colors.muted = [137 135 129]/255;   % raw measurements
colors.estimate = [42 120 214]/255; % EKF estimate
colors.faint = [195 194 183]/255;   % ground track and zoom box

fig = figure(Name="exp01: EKF track", Color="w", Position=[100 100 1100 850]);
if exist("theme", "file")
    theme(fig, "light");            % keep the saved figure light under a dark desktop theme
end
layout = tiledlayout(fig, 2, 2, TileSpacing="compact", Padding="compact");

kTurn = find(omega);
zoomCentre = X(1:3, round(mean(kTurn)))/1e3;   % middle of the turn [km]
zoomHalfWidth = 1;                             % [km]

ax = nexttile(layout);
h = plotScene(ax, radarPos, X, pMeas, xEst, colors, 8);
plotBox(ax, zoomCentre, zoomHalfWidth, colors.faint);
text(ax, radarPos(1)/1e3, radarPos(2)/1e3, radarPos(3)/1e3, "  Radar", VerticalAlignment="bottom")
text(ax, X(1, 1)/1e3, X(2, 1)/1e3, X(3, 1)/1e3, "  Start", VerticalAlignment="bottom")
view(ax, -30, 35)
title(ax, "Scenario")
lgd = legend(ax, h, Orientation="horizontal");
lgd.Layout.Tile = "north";

% Same scene, zoomed on the middle of the turn so the noise is visible.
ax = nexttile(layout);
plotScene(ax, radarPos, X, pMeas, xEst, colors, 14);
xlim(ax, zoomCentre(1) + [-1 1]*zoomHalfWidth)
ylim(ax, zoomCentre(2) + [-1 1]*zoomHalfWidth)
zlim(ax, zoomCentre(3) + [-1 1]*zoomHalfWidth)
daspect(ax, [1 1 1])
view(ax, -15, 50)
title(ax, "Zoom: middle of the turn (box at left)")

ax = nexttile(layout, [1 2]);
hold(ax, "on")
xregion(ax, t(kTurn(1)), t(kTurn(end) + 1), FaceColor=colors.muted, FaceAlpha=0.12, ...
    DisplayName="Turn");
plot(ax, t, measErr, ".", Color=colors.muted, MarkerSize=8, ...
    DisplayName=sprintf("Raw measurements (RMSE %.0f m)", rmse(measErr)));
plot(ax, t, ekfErr, "-", Color=colors.estimate, LineWidth=1.5, ...
    DisplayName=sprintf("EKF estimate (RMSE %.0f m)", rmse(ekfErr)));
hold(ax, "off")
grid(ax, "on")
xlabel(ax, "Time [s]")
ylabel(ax, "Position error [m]")
legend(ax, Location="northeast")
title(ax, "Position error")

figDir = fullfile(repoRoot, "figures");
if ~isfolder(figDir)
    mkdir(figDir);
end
exportgraphics(fig, fullfile(figDir, "exp01_ekf_track.png"), Resolution=150);

function h = plotScene(ax, radarPos, X, pMeas, xEst, colors, markerSize)
% Radar, truth, measurements and estimate on 3-D axes, in km. Returns the
% line handles in legend order.
km = 1e-3;
hold(ax, "on")
hGround = plot3(ax, km*X(1, :), km*X(2, :), zeros(1, size(X, 2)), "-", ...
    Color=colors.faint, DisplayName="Truth, ground track");
hMeas = plot3(ax, km*pMeas(1, :), km*pMeas(2, :), km*pMeas(3, :), ".", ...
    Color=colors.muted, MarkerSize=markerSize, DisplayName="Measurements");
hTruth = plot3(ax, km*X(1, :), km*X(2, :), km*X(3, :), "-", ...
    Color=colors.ink, LineWidth=1.5, DisplayName="Truth");
hEst = plot3(ax, km*xEst(1, :), km*xEst(2, :), km*xEst(3, :), "-", ...
    Color=colors.estimate, LineWidth=1.5, DisplayName="EKF estimate");
hRadar = plot3(ax, km*radarPos(1), km*radarPos(2), km*radarPos(3), "^", ...
    Color=colors.ink, MarkerFaceColor=colors.ink, MarkerSize=8, DisplayName="Radar");
hold(ax, "off")
grid(ax, "on")
xlabel(ax, "x [km]")
ylabel(ax, "y [km]")
zlabel(ax, "z [km]")
h = [hTruth, hEst, hMeas, hRadar, hGround];
end

function plotBox(ax, centre, halfWidth, color)
% Outline of the cube centre +/- halfWidth [km] that the zoom panel shows.
corners = centre + halfWidth*[-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1];
edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
hold(ax, "on")
for e = edges'
    plot3(ax, corners(1, e), corners(2, e), corners(3, e), "-", Color=color, HandleVisibility="off");
end
hold(ax, "off")
end
