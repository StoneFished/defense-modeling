function a = wrapAngle(a)
%WRAPANGLE Wrap angles to [-pi, pi).
%   a = wrapAngle(a) works elementwise. (wrapToPi does the same but needs
%   the Mapping Toolbox.)

a = mod(a + pi, 2*pi) - pi;
end
