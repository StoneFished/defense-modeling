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

% TODO: compute the target's range, azimuth and elevation relative to
% radarPos. Use atan2 so the angles land in the correct quadrant.
error("measModel:notImplemented", "Implement measModel.m.")
end
