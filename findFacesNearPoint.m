function faceIds = findFacesNearPoint(gm, P_m, radius_m, nSamples)
%FINDFACESNEARPOINT Collect geometry face IDs near a point (Monte Carlo surface sampling).
% P_m and radius_m use the same length units as gm.Vertices (meters after scaling in structuralRadiusFEA).

    if nargin < 4
        nSamples = 400;
    end

    P = P_m(:)';
    dirs = randn(nSamples, 3);
    norms = vecnorm(dirs, 2, 2);
    dirs = dirs ./ max(norms, eps);
    r = radius_m * rand(nSamples, 1).^(1/3);
    pts = P + r .* dirs;

    faceIds = zeros(nSamples, 1);
    for i = 1:nSamples
        faceIds(i) = nearestFace(gm, pts(i, 1), pts(i, 2), pts(i, 3));
    end
    faceIds = unique(faceIds);
end
