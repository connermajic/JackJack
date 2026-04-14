function faceIds = findFacesInAxisBand(gm, axis, bandLow, bandHigh, gridPerDim)
%FINDFACESINAXISBAND Collect geometry face IDs whose surface lies in an axis-aligned band.

    if nargin < 5
        gridPerDim = 6;
    end

    V = gm.Vertices;
    mn = min(V, [], 1);
    mx = max(V, [], 1);
    mn(axis) = bandLow;
    mx(axis) = bandHigh;

    xs = linspace(mn(1), mx(1), gridPerDim);
    ys = linspace(mn(2), mx(2), gridPerDim);
    zs = linspace(mn(3), mx(3), gridPerDim);
    [QX, QY, QZ] = ndgrid(xs, ys, zs);
    qpts = [QX(:), QY(:), QZ(:)];

    faceIds = zeros(size(qpts, 1), 1);
    for i = 1:size(qpts, 1)
        faceIds(i) = nearestFace(gm, qpts(i, 1), qpts(i, 2), qpts(i, 3));
    end
    faceIds = unique(faceIds);
end
