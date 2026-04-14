function stlPath = binaryVolumeToSTL(mask, p, outPath)
%BINARYVOLUMETOSTL Isosurface of binary mask -> STL file (triangle surface).

    if nargin < 3 || isempty(outPath)
        outPath = fullfile(pwd, p.tempStlName);
    end

    [nx, ny, nz] = size(mask);
    sx = p.voxelSize_mm(1);
    sy = p.voxelSize_mm(2);
    sz = p.voxelSize_mm(3);
    ox = p.volumeOrigin_mm(1);
    oy = p.volumeOrigin_mm(2);
    oz = p.volumeOrigin_mm(3);

    xg = ox + (0:nx - 1) * sx;
    yg = oy + (0:ny - 1) * sy;
    zg = oz + (0:nz - 1) * sz;

    [X, Y, Z] = ndgrid(xg, yg, zg);
    B = double(mask);

    fv = isosurface(X, Y, Z, B, p.isovalue);
    if isempty(fv.vertices)
        error('binaryVolumeToSTL:Empty', 'Isosurface produced no vertices. Check segmentation and isovalue.');
    end

    try
        fv = reducepatch(fv, 1.0);
    catch
        % reducepatch optional
    end

    TR = triangulation(fv.faces, fv.vertices);
    stlwrite(TR, outPath);
    stlPath = outPath;
end
