function H = measJacobian(x, radarPos)
%MEASJACOBIAN Jacobian of the radar measurement function, H = dh/dx (3x6).
%   H = measJacobian(x, radarPos) linearizes measModel at the state
%   x = [px; py; pz; vx; vy; vz]. Row i is the gradient of the i-th
%   measurement (range, azimuth, elevation) with respect to x.
%
%   Derive it by hand. tests/MeasurementModelTest.m compares it against
%   central finite differences of measModel. H is undefined when the target
%   is directly above or below the radar (zero horizontal range); you do not
%   need to handle that case.
%
%   See also measModel, ekf.

arguments
    x (6,1) double
    radarPos (3,1) double
end

% TODO: fill in H analytically.
error("measJacobian:notImplemented", "Implement measJacobian.m.")
end
