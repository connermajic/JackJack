function A = faceAreasFromGeometry(gm)
%FACEAREASFROMGEOMETRY Approximate each geometry face area via surface triangulation.

    try
        TR = triangulation(gm);
    catch ME
        error('faceAreasFromGeometry:Triangulation', ...
            'Could not build triangulation from geometry: %s', ME.message);
    end

    P = TR.Points;
    C = TR.ConnectivityList;
    nf = gm.NumFaces;
    A = zeros(nf, 1);

    for i = 1:size(C, 1)
        triPts = P(C(i, :), :);
        c = mean(triPts, 1);
        fid = nearestFace(gm, c(1), c(2), c(3));
        v1 = triPts(2, :) - triPts(1, :);
        v2 = triPts(3, :) - triPts(1, :);
        triArea = 0.5 * norm(cross(v1, v2));
        if fid >= 1 && fid <= nf
            A(fid) = A(fid) + triArea;
        end
    end
end
