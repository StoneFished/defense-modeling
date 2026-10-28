function e2 = nees(X, xEst, PEst)
%NEES Normalized estimation error squared of a filter run.
%   e2 = nees(X, xEst, PEst) returns, for each column k, the squared error
%   of the estimate xEst(:,k) against the truth X(:,k), normalized by the
%   filter's covariance PEst(:,:,k):
%
%       e2(k) = e'*inv(PEst(:,:,k))*e,   e = xEst(:,k) - X(:,k).
%
%   If the filter's covariance is honest, the average of e2(k) over many
%   Monte Carlo runs equals the state dimension (6 for [p; v]). Much
%   larger means the filter is overconfident; smaller, conservative.
%
%   See also ekf, immFilter.

arguments
    X (:,:) double
    xEst (:,:) double {mustBeEqualSize(X, xEst)}
    PEst (:,:,:) double
end

e = xEst - X;
e2 = zeros(1, size(X, 2));
for k = 1:size(X, 2)
    e2(k) = e(:, k)' / PEst(:, :, k) * e(:, k);
end
end

function mustBeEqualSize(a, b)
if ~isequal(size(a), size(b))
    error("nees:sizeMismatch", "The truth and the estimates must be the same size.")
end
end
