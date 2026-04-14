function [V, meta] = loadConvertedVolume(imagePath, opts)
%LOADCONVERTEDVOLUME Load a converted scan volume (3-D array) from disk.
%
% Supported:
%   *.mat          — loads a 3-D variable (default: largest 3-D numeric array)
%   *.nii, *.nii.gz — uses niftiread (Image Processing Toolbox / Medical Imaging)
%   *.tif, *.tiff  — multi-page TIFF stack (reads all frames into 3rd dim)
%   *.isq          — raw ISQ volume payload (with optional header parsing)
%
% opts (optional struct):
%   .matVariableName — char, force variable name for .mat files
%   .isqSize        — [nx ny nz], required if not inferable from ISQ header text
%   .isqDataType    — fread type (default 'int16')
%   .isqHeaderBytes — byte offset to payload (default 512)
%   .isqByteOrder   — fopen machine format, e.g. 'ieee-le' (default) or 'ieee-be'

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
        case '.isq'
            V = loadRawIsqVolume(imagePath, opts);
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

function V = loadRawIsqVolume(imagePath, opts)
    if isfield(opts, 'isqDataType') && ~isempty(opts.isqDataType)
        dataType = opts.isqDataType;
    else
        dataType = 'int16';
    end

    if isfield(opts, 'isqHeaderBytes') && ~isempty(opts.isqHeaderBytes)
        headerBytes = opts.isqHeaderBytes;
    else
        headerBytes = 512;
    end

    if isfield(opts, 'isqByteOrder') && ~isempty(opts.isqByteOrder)
        byteOrder = opts.isqByteOrder;
    else
        byteOrder = 'ieee-le';
    end

    fileInfo = dir(imagePath);
    if isempty(fileInfo)
        error('loadConvertedVolume:IsqNotFound', 'Unable to stat ISQ file: %s', imagePath);
    end

    if isfield(opts, 'isqSize') && ~isempty(opts.isqSize)
        volSize = opts.isqSize;
    else
        volSize = inferIsqSizeFromHeader(imagePath);
        if isempty(volSize)
            error('loadConvertedVolume:IsqSizeRequired', ...
                ['Could not infer ISQ volume size. Provide opts.isqSize = [nx ny nz] ', ...
                 '(via paramOverride.isqSize in radiusFEA_run).']);
        end
    end

    if ~isnumeric(volSize) || numel(volSize) ~= 3 || any(volSize <= 0)
        error('loadConvertedVolume:IsqBadSize', 'opts.isqSize must be [nx ny nz] with positive values.');
    end
    volSize = double(volSize(:)');

    fid = fopen(imagePath, 'r', byteOrder);
    if fid < 0
        error('loadConvertedVolume:IsqOpen', 'Unable to open ISQ file: %s', imagePath);
    end
    cleaner = onCleanup(@() fclose(fid));

    status = fseek(fid, headerBytes, 'bof');
    if status ~= 0
        error('loadConvertedVolume:IsqSeek', 'Unable to seek to header offset %d in %s', headerBytes, imagePath);
    end

    data = fread(fid, prod(volSize), ['*' dataType]);
    if numel(data) ~= prod(volSize)
        error('loadConvertedVolume:IsqRead', ...
            'Expected %d voxels from ISQ payload, read %d.', prod(volSize), numel(data));
    end

    expectedBytes = headerBytes + prod(volSize) * bytesPerSample(dataType);
    if expectedBytes > fileInfo.bytes
        error('loadConvertedVolume:IsqSizeMismatch', ...
            'Header/size/type imply %d bytes, but file has %d bytes.', expectedBytes, fileInfo.bytes);
    end

    V = reshape(data, volSize);
end

function n = bytesPerSample(dataType)
    dt = lower(strtrim(dataType));
    switch dt
        case {'int8', 'uint8', 'char', 'schar', 'uchar'}
            n = 1;
        case {'int16', 'uint16', 'short', 'ushort'}
            n = 2;
        case {'int32', 'uint32', 'int', 'uint', 'single', 'float32'}
            n = 4;
        case {'int64', 'uint64', 'double', 'float64'}
            n = 8;
        otherwise
            error('loadConvertedVolume:IsqType', 'Unsupported opts.isqDataType: %s', dataType);
    end
end

function volSize = inferIsqSizeFromHeader(imagePath)
    volSize = [];

    fid = fopen(imagePath, 'r', 'ieee-le');
    if fid < 0
        return;
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    headerChars = fread(fid, 4096, '*char')';
    if isempty(headerChars)
        return;
    end

    tokX = regexp(headerChars, 'dimx_p\s*[:=]\s*(\d+)', 'tokens', 'once');
    tokY = regexp(headerChars, 'dimy_p\s*[:=]\s*(\d+)', 'tokens', 'once');
    tokZ = regexp(headerChars, 'dimz_p\s*[:=]\s*(\d+)', 'tokens', 'once');
    if ~isempty(tokX) && ~isempty(tokY) && ~isempty(tokZ)
        volSize = [str2double(tokX{1}), str2double(tokY{1}), str2double(tokZ{1})];
        return;
    end

    tokX = regexp(headerChars, 'dimx\s*[:=]\s*(\d+)', 'tokens', 'once');
    tokY = regexp(headerChars, 'dimy\s*[:=]\s*(\d+)', 'tokens', 'once');
    tokZ = regexp(headerChars, 'dimz\s*[:=]\s*(\d+)', 'tokens', 'once');
    if ~isempty(tokX) && ~isempty(tokY) && ~isempty(tokZ)
        volSize = [str2double(tokX{1}), str2double(tokY{1}), str2double(tokZ{1})];
    end
end
