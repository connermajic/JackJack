function [V, meta] = loadConvertedVolume(imagePath, opts)
%LOADCONVERTEDVOLUME Load a converted scan volume (3-D array) from disk.
%
% Supported:
%   *.mat          — loads a 3-D variable (default: largest 3-D numeric array)
%   *.nii, *.nii.gz — uses niftiread (Image Processing Toolbox / Medical Imaging)
%   *.tif, *.tiff  — multi-page TIFF stack (reads all frames into 3rd dim)
%
% opts (optional struct):
%   .matVariableName — char, force variable name for .mat files

    if nargin < 2
        opts = struct();
    end

    meta = struct('path', imagePath, 'class', '', 'size', []);

    if ~isfile(imagePath)
        error('loadConvertedVolume:NotFound', 'File not found: %s', imagePath);
    end

    [~, ~, ext] = fileparts(imagePath);
    ext = lower(ext);
    lowerPath = lower(imagePath);

    if endsWith(lowerPath, '.nii.gz')
        if exist('niftiread', 'file') ~= 2
            error('loadConvertedVolume:MissingNifti', 'niftiread not found. Use a .mat volume or install NIfTI support.');
        end
        V = niftiread(imagePath);
        meta.path = imagePath;
        meta.class = class(V);
        V = double(V);
        if ndims(V) ~= 3
            error('loadConvertedVolume:Dims', 'Expected a 3-D volume, got size %s.', mat2str(size(V)));
        end
        meta.size = size(V);
        return;
    end

    switch ext
        case '.mat'
            S = load(imagePath);
            if isfield(opts, 'matVariableName') && ~isempty(opts.matVariableName)
                vn = opts.matVariableName;
                if ~isfield(S, vn)
                    error('loadConvertedVolume:BadMat', 'Variable "%s" not in %s', vn, imagePath);
                end
                V = S.(vn);
            else
                V = pickLargest3dNumeric(S);
            end
        case '.nii'
            if exist('niftiread', 'file') ~= 2
                error('loadConvertedVolume:MissingNifti', 'niftiread not found.');
            end
            V = niftiread(imagePath);
        case {'.tif', '.tiff'}
            info = imfinfo(imagePath);
            n = numel(info);
            if n == 1
                V = imread(imagePath);
                if size(V, 3) > 1
                    % RGB etc. — take first channel
                    V = V(:, :, 1);
                end
            else
                V = zeros(info(1).Height, info(1).Width, n, 'like', imread(imagePath, 1));
                for k = 1:n
                    Vk = imread(imagePath, k);
                    if size(Vk, 3) > 1
                        Vk = Vk(:, :, 1);
                    end
                    V(:, :, k) = Vk;
                end
            end
        otherwise
            error('loadConvertedVolume:Unsupported', 'Unsupported extension: %s', ext);
    end

    if ~isnumeric(V) && ~islogical(V)
        error('loadConvertedVolume:Type', 'Volume must be numeric or logical.');
    end

    V = double(V);
    if ndims(V) ~= 3
        error('loadConvertedVolume:Dims', 'Expected a 3-D volume, got size %s.', mat2str(size(V)));
    end

    meta.class = class(V);
    meta.size = size(V);
end

function V = pickLargest3dNumeric(S)
    names = fieldnames(S);
    best = [];
    bestName = '';
    for i = 1:numel(names)
        A = S.(names{i});
        if isnumeric(A) && ndims(A) == 3
            n = numel(A);
            if isempty(best) || n > best
                best = n;
                bestName = names{i};
            end
        end
    end
    if isempty(bestName)
        error('loadConvertedVolume:No3D', 'No 3-D numeric array found in MAT file.');
    end
    V = S.(bestName);
end
