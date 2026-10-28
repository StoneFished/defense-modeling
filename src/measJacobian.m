function H = measJacobian(x, radarPos)
%MEASJACOBIAN Jacobian of the radar measurement function, H = dh/dx (3x6).
%   H = measJacobian(x, radarPos) linearizes measModel at the state
%   x = [px; py; pz; vx; vy; vz]. Row i is the gradient of the i-th
%   measurement (range, azimuth, elevation) with respect to x.
%
%   H is derived by hand; tests/MeasurementModelTest.m compares it against
%   central finite differences of measModel. H is undefined when the target
%   is directly above or below the radar (zero horizontal range); that case
%   is not handled.
%
%   See also measModel, ekf.

arguments
    x (6,1) double
    radarPos (3,1) double
end

d = x(1:3) - radarPos;              % radar-to-target vector
rho2 = d(1)^2 + d(2)^2;             % horizontal range, squared
rho = sqrt(rho2);
r2 = rho2 + d(3)^2;                 % range, squared
r = sqrt(r2);

% Rows: range = norm(d), azimuth = atan2(dy, dx), elevation = atan2(dz, rho).
% The measurements do not depend on velocity, so columns 4:6 stay zero.
H = zeros(3, 6);
H(1, 1:3) = d'/r;
H(2, 1:2) = [-d(2), d(1)]/rho2;
H(3, 1:3) = [-d(1)*d(3), -d(2)*d(3), rho2]/(r2*rho);
end
