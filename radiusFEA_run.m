function result = radiusFEA_run(imagePath, paramOverride)
%RADIUSFEA_RUN End-to-end homogeneous distal-radius FEA from a converted volume.
%
% Usage:
%   result = radiusFEA_run('C:\data\radius.mat');
%   result = radiusFEA_run('scan.nii.gz', struct('forceTotal_N',[0 0 -800], 'loadOffset_mm',[3 0 0]));
%
% Key outputs (see defaultRadiusFEAParams):
%   result.loadPointTarget_mm — target point for the off-axis load patch (image/world axes, mm)
%   result.surfaceTraction_Pa — uniform traction on the patch (vector) giving the requested total force
%   result.solution — PDE solution (displacement, stress)
%
% Requires Partial Differential Equation Toolbox. Image segmentation uses Image Processing Toolbox
% helpers when available (graythresh, imfill, imopen).

    if nargin < 1 || isempty(imagePath)
        error('radiusFEA_run:Input', 'Provide the path to a converted volume (.mat, .nii, .nii.gz, .tif, .isq).');
    end

    p = defaultRadiusFEAParams();
    if nargin >= 2 && ~isempty(paramOverride)
        p = mergeStruct(p, paramOverride);
    end

    opts = struct();
    if nargin >= 2 && isstruct(paramOverride) && isfield(paramOverride, 'matVariableName')
        opts.matVariableName = paramOverride.matVariableName;
    end
    if nargin >= 2 && isstruct(paramOverride) && isfield(paramOverride, 'isqSize')
        opts.isqSize = paramOverride.isqSize;
    end
    if nargin >= 2 && isstruct(paramOverride) && isfield(paramOverride, 'isqDataType')
        opts.isqDataType = paramOverride.isqDataType;
    end
    if nargin >= 2 && isstruct(paramOverride) && isfield(paramOverride, 'isqHeaderBytes')
        opts.isqHeaderBytes = paramOverride.isqHeaderBytes;
    end
    if nargin >= 2 && isstruct(paramOverride) && isfield(paramOverride, 'isqByteOrder')
        opts.isqByteOrder = paramOverride.isqByteOrder;
    end

    [V, meta] = loadConvertedVolume(imagePath, opts);
    fprintf('Loaded volume size %s\n', mat2str(meta.size));

    [mask, segInfo] = segmentBoneVolume(V, p);
    fprintf('Segmentation threshold = %.6g (Otsu=%d)\n', segInfo.threshold, segInfo.usedOtsu);

    stlPath = binaryVolumeToSTL(mask, p, fullfile(pwd, p.tempStlName));
    fprintf('Wrote STL: %s\n', stlPath);

    result = structuralRadiusFEA(stlPath, p);
    result.volumeMeta = meta;
    result.segmentation = segInfo;

    fprintf('Load application target (mm, same axes as volume): [%.4f, %.4f, %.4f]\n', ...
        result.loadPointTarget_mm(1), result.loadPointTarget_mm(2), result.loadPointTarget_mm(3));
    fprintf('Total force (N): [%.4f, %.4f, %.4f]\n', result.forceTotal_N(1), result.forceTotal_N(2), result.forceTotal_N(3));
    fprintf('Load patch area (m^2): %.6e\n', result.loadPatchArea_m2);

    if p.writeResultMat
        save(p.resultMatFile, 'result', 'p', '-v7.3');
        fprintf('Saved %s\n', p.resultMatFile);
    end
end

function out = mergeStruct(a, b)
    out = a;
    if isempty(b) || ~isstruct(b)
        return;
    end
    fn = fieldnames(b);
    for i = 1:numel(fn)
        out.(fn{i}) = b.(fn{i});
    end
end
