function z = measModel(x, radarPos)
%MEASMODEL Radar measurement function h(x): range, azimuth, elevation.
%   z = measModel(x, radarPos) maps the state x = [px; py; pz; vx; vy; vz]
%   to z = [range; azimuth; elevation] as seen from a radar at radarPos
%   (3x1, same Cartesian frame as x).
%
%   Conventions (tests/MeasurementModelTest.m checks these):
%       range      distance from the radar to the target [m]
%       azimuth    angle in the x-y plane from +x toward +y, in [-pi, pi]
%                  [rad] (a math angle, not a compass bearing)
%       elevation  angle above the x-y plane, in [-pi/2, pi/2] [rad]
%
%   See also measJacobian, meas2cart, generateMeasurements.

arguments
    x (6,1) double
    radarPos (3,1) double
end

d = x(1:3) - radarPos;              % radar-to-target vector
z = [norm(d); atan2(d(2), d(1)); atan2(d(3), hypot(d(1), d(2)))];
end
