function p = meas2cart(z, radarPos)
%MEAS2CART Convert radar measurements to Cartesian positions.
%   p = meas2cart(z, radarPos) converts the columns of z (3 x N, rows
%   [range; azimuth; elevation], conventions as in measModel) taken by a
%   radar at radarPos into positions p (3 x N) in the common Cartesian frame.
%
%   See also measModel.

arguments
    z (3,:) double
    radarPos (3,1) double
end

r = z(1, :);
az = z(2, :);
el = z(3, :);
p = radarPos + [r.*cos(el).*cos(az); r.*cos(el).*sin(az); r.*sin(el)];
end
