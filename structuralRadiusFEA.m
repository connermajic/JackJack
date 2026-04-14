function result = structuralRadiusFEA(stlPath, p)
%STRUCTURALRADIUSFEA Homogeneous linear elastic FEA with off-axis distal load patch.

    if exist('createpde', 'file') ~= 2
        error('structuralRadiusFEA:NoPDE', 'Partial Differential Equation Toolbox is required (createpde).');
    end

    result = struct();
    result.stlPath = stlPath;

    model = createpde('structural', 'static-solid');

    try
        importGeometry(model, stlPath, 'FeatureAngle', p.stlFeatureAngle);
    catch ME1
        try
            importGeometry(model, stlPath);
        catch ME2
            error('structuralRadiusFEA:Import', 'importGeometry failed: %s / %s', ME1.message, ME2.message);
        end
    end

    gm = model.Geometry;

    % mm -> m for SI (N, Pa)
    try
        scale(model.Geometry, 1e-3);
    catch
        try
            scale(model.Geometry, [1e-3, 1e-3, 1e-3]);
        catch ME
            error('structuralRadiusFEA:Scale', 'Could not scale geometry to meters: %s', ME.message);
        end
    end
    gm = model.Geometry;

    axis = p.axisForLength;
    V = gm.Vertices;
    coords = V(:, axis);
    cmin = min(coords);
    cmax = max(coords);
    extent = cmax - cmin;

    if isempty(p.proximalBand_mm)
        bandP = 0.03 * extent;
    else
        bandP = p.proximalBand_mm * 1e-3;
    end
    if isempty(p.distalBand_mm)
        bandD = 0.03 * extent;
    else
        bandD = p.distalBand_mm * 1e-3;
    end

    faceProx = findFacesInAxisBand(gm, axis, cmin, cmin + bandP);
    faceDist = findFacesInAxisBand(gm, axis, cmax - bandD, cmax);

    if isempty(faceProx)
        error('structuralRadiusFEA:NoProximal', 'No proximal faces found; check axisForLength or bands.');
    end
    if isempty(faceDist)
        error('structuralRadiusFEA:NoDistal', 'No distal faces found; check axisForLength or bands.');
    end

    % Distal centroid (meters) from vertices near distal end
    axv = V(:, axis);
    distalVertMask = axv >= (cmax - bandD);
    distalCentroid = mean(V(distalVertMask, :), 1);

    F = p.forceTotal_N(:);
    off = projectOffsetOrthogonalToForce(p.loadOffset_mm * 1e-3, F);
    loadPoint_m = distalCentroid + off;

    faceLoad = findFacesNearPoint(gm, loadPoint_m, p.loadPatchRadius_mm * 1e-3, 500);
    faceLoad = intersect(faceLoad, faceDist);
    if isempty(faceLoad)
        faceLoad = faceDist;
        warning('structuralRadiusFEA:ExpandedLoad', ...
            'No faces in offset patch; applying uniform traction on full distal band.');
    end

    A = faceAreasFromGeometry(gm);
    Aload = sum(A(faceLoad));
    if Aload < eps
        error('structuralRadiusFEA:ZeroArea', 'Zero load-patch area; segmentation or bands may be wrong.');
    end

    traction_Pa = F / Aload;

    structuralProperties(model, 'YoungsModulus', p.youngModulus_Pa, 'PoissonsRatio', p.poissonRatio);
    structuralBC(model, 'Face', faceProx, 'Constraint', 'fixed');
    structuralBoundaryLoad(model, 'Face', faceLoad, 'SurfaceTraction', traction_Pa);

    if isempty(p.meshHmax_mm)
        hmax = 0.15 * min(max(V) - min(V));
    else
        hmax = p.meshHmax_mm * 1e-3;
    end

    if isempty(p.meshHmin_mm)
        generateMesh(model, 'Hmax', hmax, 'Hgrad', p.meshHgrad);
    else
        generateMesh(model, 'Hmax', hmax, 'Hmin', p.meshHmin_mm * 1e-3, 'Hgrad', p.meshHgrad);
    end

    R = solve(model);

    result.model = model;
    result.solution = R;
    result.faceIdsProximal = faceProx;
    result.faceIdsDistalBand = faceDist;
    result.faceIdsLoad = faceLoad;
    result.loadPatchArea_m2 = Aload;
    result.surfaceTraction_Pa = traction_Pa;
    result.distalCentroid_m = distalCentroid;
    result.loadPointTarget_m = loadPoint_m;
    result.loadPointTarget_mm = loadPoint_m * 1e3;
    result.forceTotal_N = F;
    result.meshStats = struct('NumNodes', model.Mesh.NumNodes, 'NumElements', model.Mesh.NumElements);

    if p.showPlots
        try
            figure('Name', 'Displacement magnitude (m)');
            pdeplot3D(model.Mesh, 'ColorMapData', R.Displacement.Magnitude);
            title('Displacement magnitude');
            colorbar;
        catch
        end
        try
            figure('Name', 'von Mises stress (Pa)');
            if isfield(R, 'VonMisesStress')
                pdeplot3D(model.Mesh, 'ColorMapData', R.VonMisesStress);
            elseif isfield(R, 'vonMisesStress')
                pdeplot3D(model.Mesh, 'ColorMapData', R.vonMisesStress);
            end
            title('von Mises stress');
            colorbar;
        catch
        end
    end
end
